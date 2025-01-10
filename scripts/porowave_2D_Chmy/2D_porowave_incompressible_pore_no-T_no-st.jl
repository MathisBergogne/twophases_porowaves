using Chmy
using KernelAbstractions
using Printf
using CairoMakie

@views avx(A) = 0.5 .* (A[1:end-1,:] .+ A[2:end,:])
@views avy(A) = 0.5 .* (A[:,1:end-1] .+ A[:,2:end])

@views function Test_compute_qD()
	lx   = 2500      # longeur du modele en x
	ly   = 5000     # longeur du modele en y
	k_r  = 1e-11     # permeabilite initial m2 / viscosité fluide Pa s
	η_f	 = 1e4       # viscosité fluide Pa s
	η_r  = 1e17      # matrix bulk viscosity
	ϕ_bg = 0.01      # porosité du fond
	loc_ϕ = 0.85	 # localisation of the intitial porosity on the box (0-top/1-bottom)
	w	 = 500		 # shape aspect of the initial melt 
	npow = 3         # exposant de la porosité de l'équation Carman–Kozeny
	ϕA   = 0.1       # différence de porosité entre le fond et une bulle de porosité
	ρs		= 2700							# matrix volumetric mass (kg.m-3)
	ρf		= 2300							# mvolumetric mass fluid (kg.m-3)
	g		= 9.81							# gravity (m.s-2)
	Bf   = 1e-1       # compressibilité du fluide 
	# independent physics
	lc		= sqrt(k_r / η_f * η_r / ϕ_bg)	# m
	Δρg		= (ρs - ρf) * g					# Pa.m-1
	η_ϕbg	= η_r / ϕ_bg
	psc		= Δρg * lc				# Pa
	tsc		= η_ϕbg / psc	# s
	println("lc : ",lc,";	tsc : ",tsc)
	# numerics
	nx   = 100       # nombre de noeuds de modèles en x
	ny   = 200       # nombre de noeuds de modèles en y
	dx   = lx / nx   # pas d'espace en x
	dy   = ly / ny   # pas d'espace en y
	xc   = LinRange((-lx+dx)/2,(lx-dx)/2,nx)   # coordonnées en x des noeuds du modèles
	yc   = LinRange(-ly+(dy)/2,-dy/2,ny)   # coordonnées en y des noeuds du modèles
	x_tick = [-10^round(log10(xc[end])),0,10^round(log10(xc[end]))]
	nt   = 1e7      # nombre de pas de temps
	dt = 1e-3 * tsc
	nvis = nt / 100
	year2sec = 365*24*3600
	### initialisation ###
	backend=CPU()
	arch = Arch(backend)
	nxy=(nx, ny)
    # geometry
    grid   = UniformGrid(arch; origin=(-lx/2, -ly), extent=(lx, ly), dims=nxy)
    launch = Launcher(arch, grid; outer_width=(16, 8))
	# fields
	ϕ = Field(backend, grid, Center())
	Pe = Field(backend, grid, Center())
	k_rg = Field(backend, grid, Center())
	η_fg = Field(backend, grid, Center())
	η_rg = Field(backend, grid, Center())
	dϕdt =  Field(backend, grid, Center())
	dPedt =  Field(backend, grid, Center())
    qD = VectorField(backend, grid)
    # initial conditions
	set!(ϕ, grid, (x, y) -> ϕ_bg + ϕA * exp(-((x * loc_ϕ) / w) ^ 2 -((y + ly * loc_ϕ) / w) ^ 2))
	bc!(arch, grid, ϕ => Neumann(); exchange=ϕ)
    set!(Pe, grid, (_, _) -> 0.0)
	bc!(arch, grid, Pe => Neumann(); exchange=Pe)
    set!(dPedt, grid, (_, _) -> 0.0)
    set!(dϕdt, grid, (_, _) -> 0.0)
	set!(k_rg, grid, (_, _) -> k_r)
    set!(η_fg, grid, (_, _) -> η_f)
    set!(η_rg, grid, (_, _) -> η_r)
	# visualisation
	ploting(xc,yc,ϕ,qDx,qDy,dϕdt,Pe,dPedt,x_tick,real_time)
	#loop
	for it = 1:nt 
		# compute qD.x and qD.y
		launch(arch, grid, compute_qD! => (qD, ϕ, Pe, k_r, η_fg, Δρg, npow, grid))
	    launch(arch, grid, update_pressure! => (qD, ϕ, Pe, dϕdt, dPedt, Bf, dt, grid); bc=batch(grid, Pe => Neumann(); exchange=Pe))
	    launch(arch, grid, update_porosity! => (ϕ, Pe, dϕdt, η_rg, dt, grid); bc=batch(grid, ϕ => Neumann(); exchange=ϕ))
		if (it % (nvis/10)) == 0
			println("iter :",it)
		end
			# visualisation
	    if (it % nvis) == 0
			print("MAX : \n dPedt : ",maximum(dPedt),"; ϕ : ",maximum(ϕ),"\n")
			print("MIN : \n dPedt : ",minimum(dPedt),"; ϕ : ",maximum(ϕ),"\n")
			print("Bf : ", Bf,"\n")
			# update visualisation
	    	KernelAbstractions.synchronize(backend)
			ploting(xc,yc,ϕ,qDx,qDy,dϕdt,Pe,dPedt,x_tick,real_time)
		end
	end
end

@kernel inbounds = true function compute_qD!(qD, ϕ, Pe, k_r,η_fg, Δρg, npow, g::StructuredGrid, O)
    I = @index(Global, NTuple)
    I = I + O
	# qD.x = (k_r/η_fg) * ϕ^npow * ( ∂x(Pe, g, I...) )
	qD.x[I...] = (k_r/lerp(η_fg,location(qD.x), g, I...)) * lerp(ϕ,location(qD.x), g, I...).^npow * ( ∂x(Pe, g, I...) )
    # qD.y = (k_r/η_fg) * ϕ^npow * ( ∂y(Pe, g, I...) - Δρg)
	qD.y[I...] = (k_r/lerp(η_fg,location(qD.y), g, I...)) * lerp(ϕ,location(qD.y), g, I...).^npow * ( ∂y(Pe, g, I...) - Δρg)
end

@kernel inbounds = true function update_porosity!(ϕ, Pe, dϕdt, η_rg, dt, g::StructuredGrid, O)
    I = @index(Global, NTuple)
    I = I + O
	dϕdt[I...] = -(1 - ϕ[I...])* (ϕ[I...] / η_rg[I...]) * Pe[I...]
	ϕ[I...] += dt * dϕdt[I...]
end

@kernel inbounds = true function update_pressure!(qD, ϕ, Pe, dϕdt, dPedt, Bf, dt, g::StructuredGrid, O)
    I = @index(Global, NTuple)
    I = I + O
	dPedt[I...] = (1 / (ϕ[I...] .* Bf)) * ( 1 / (1 .- ϕ[I...]) * dϕdt[I...] - (∂x(qD.x, g, I...) + (∂y(qD.y, g, I...))))
	Pe[I...] += dt * dPedt[I...]
end

@views function def_ticks(A)
	M = (max(abs(minimum(A)),maximum(A)))
	if M != 0 
		P = round(log10(M),RoundDown)
		L = round(M/10^P,RoundDown)*10^P
		if isnan(L) 
			print("M : ",M,"P : ",P,"L :",L,"\n")
		end
		A_ticks = -L:L/2:L 
	else
		P = 0
		L = 0 
		A_ticks = -0.5:0.5:0.5 
	end
	return A_ticks
end

@views function ploting(xc,yc,ϕ,qDx,qDy,dϕdt,Pe,dPedt,x_tick,real_time)
	fig = Figure(fontsize = 11; size=(600, 600))
	axs = (ϕ     = Axis(fig[1, 1][1, 1]; aspect=DataAspect(), xlabel="x", title="ϕ",    xticks = x_tick, xtickformat = "{:.0f}", xticklabelrotation = pi/4, ylabel="Depth (m)"),
		   qDx   = Axis(fig[1, 2][1, 1]; aspect=DataAspect(), xlabel="x", title="qDx",  xticks = x_tick, xtickformat = "{:.0f}", xticklabelrotation = pi/4),
		   qDy   = Axis(fig[1, 3][1, 1]; aspect=DataAspect(), xlabel="x", title="qDy",  xticks = x_tick, xtickformat = "{:.0f}", xticklabelrotation = pi/4),
		   dϕdt  = Axis(fig[2, 1][1, 1]; aspect=DataAspect(), xlabel="x", title="dϕdt", xticks = x_tick, xtickformat = "{:.0f}", xticklabelrotation = pi/4, ylabel="Depth (m)"),
		   Pe    = Axis(fig[2, 2][1, 1]; aspect=DataAspect(), xlabel="x", title="Pe",   xticks = x_tick, xtickformat = "{:.0f}", xticklabelrotation = pi/4),
		   dPedt = Axis(fig[2, 3][1, 1]; aspect=DataAspect(), xlabel="x", title="dPedt",xticks = x_tick, xtickformat = "{:.0f}", xticklabelrotation = pi/4))
	plt = (ϕ     = heatmap!(axs.ϕ  ,   xc,yc, ϕ     |> Array; colormap=:turbo, colorrange = (0.0, max(0.15,maximum(ϕ)))),
		   qDx   = heatmap!(axs.qDx,   xc,yc, qDx   |> Array; colormap=:turbo),
		   qDy   = heatmap!(axs.qDy,   xc,yc, qDy   |> Array; colormap=:turbo), 
		   dϕdt  = heatmap!(axs.dϕdt,  xc,yc, dϕdt  |> Array; colormap=:turbo),
		   Pe    = heatmap!(axs.Pe,    xc,yc, Pe    |> Array; colormap=:turbo),
		   dPedt = heatmap!(axs.dPedt, xc,yc, dPedt |> Array; colormap=:turbo))
	Colorbar(fig[1, 1][1, 2], plt.ϕ,     ticks = 0.0:0.05:1)
	Colorbar(fig[1, 2][1, 2], plt.qDx,   ticks = def_ticks(qDx), tickformat = "{:.0e}")
	Colorbar(fig[1, 3][1, 2], plt.qDy,   ticks = def_ticks(qDy), tickformat = "{:.0e}")
	Colorbar(fig[2, 1][1, 2], plt.dϕdt,  ticks = def_ticks(dϕdt), tickformat = "{:.0e}")
	Colorbar(fig[2, 2][1, 2], plt.Pe,    ticks = def_ticks(Pe), tickformat = "{:.0e}")
	Colorbar(fig[2, 3][1, 2], plt.dPedt, ticks = def_ticks(dPedt), tickformat = "{:.0e}")
	Label(fig[0, :], text = @sprintf("time : %1.1e y",real_time), fontsize = 14)
   	display(fig)
end

Test_compute_qD()
