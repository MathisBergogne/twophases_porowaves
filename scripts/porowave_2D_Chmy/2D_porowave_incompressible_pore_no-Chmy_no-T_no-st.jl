using CairoMakie, Printf

@views avx(A) = 0.5 .* (A[1:end-1,:] .+ A[2:end,:])
@views avy(A) = 0.5 .* (A[:,1:end-1] .+ A[:,2:end])

@views function porowave_1D_pore(do_visu)
	lx   = 2500       # longeur du modele en x
	ly   = 5000      # longeur du modele en y
	k_r  = 1e-11     # permeabilite initial m2 / viscosité fluide Pa s
	η_f	 = 1e4       # viscosité fluide Pa s
	η_r  = 1e17      # matrix bulk viscosity
	ϕ_bg = 0.01      # porosité du fond
	loc_ϕ = 0.85	 # localisation of the intitial porosity on the box (0-top/1-bottom)
	w	 = 500		 # shape aspect of the initial melt 
	npow = 3         # exposant de la porosité de l'équation Carman–Kozeny
	ϕA   = 0.1       # différence de porosité entre le fond et une bulle de porosité
	ρs	 = 2700							# matrix volumetric mass (kg.m-3)
	ρf	 = 2300							# mvolumetric mass fluid (kg.m-3)
	g	 = 9.81							# gravity (m.s-2)
	Bf   = 1e-1      # compressibilité du fluide 
	# independent physics
	lc		= sqrt(k_r / η_f * η_r / ϕ_bg)	# m
	Δρg		= (ρs - ρf) * g					# Pa.m-1
	η_ϕbg	= η_r / ϕ_bg
	psc		= Δρg * lc				# Pa
	tsc		= η_ϕbg / psc	# s
	println("lc : ",lc,";	tsc : ",tsc)
	# numerics
	nx   = 100       # nombre de noeuds de modèles en x
	ny   = 200      # nombre de noeuds de modèles en y
	dx   = lx / nx   # pas d'espace en x
	dy   = ly / ny   # pas d'espace en y
	xc   = LinRange((-lx+dx)/2,(lx-dx)/2,nx)   # coordonnées en x des noeuds du modèles
	yc   = LinRange(-ly+(dy)/2,-dy/2,ny)   # coordonnées en y des noeuds du modèles
	x_tick = [-10^round(log10(xc[end])),0,10^round(log10(xc[end]))]
	year2sec = 365*24*3600
	nt   = 1e7      # nombre de pas de temps
	dt = 1e-3 * tsc 
	nvis = nt / 100
	### initialisation ###
	# mecanics
	ϕ		= @. ϕ_bg + ϕA * exp(-((xc * loc_ϕ) / w) ^ 2 -((yc' + ly * loc_ϕ) / w) ^ 2)
	Pe     = zeros(nx,ny)
	qDx    = zeros(nx - 1,ny - 2)
	qDy    = zeros(nx - 2,ny - 1)
	dPedt  = zeros(nx,ny)
	dϕdt   = zeros(nx,ny)
	if do_visu ploting_ϕ(xc,yc,ϕ,qDx,qDy,dϕdt,Pe,dPedt,x_tick,0.0) end
	#loop
	for it = 1:nt
	    compute_meca(qDx, qDy,dϕdt, dPedt, Pe, ϕ, k_r, η_f, η_r, Bf, npow, dx, dy, dt, Δρg)
		if (it % (nvis/10)) == 0
			println("iter :",it)
		end
	    # visualisation
	    if (it % nvis) == 0 && do_visu || isnan(ϕ[2])
			ploting_ϕ(xc,yc,ϕ,qDx,qDy,dϕdt,Pe,dPedt,x_tick,dt*it/year2sec)
	    end
		if isnan(ϕ[2]); 
			println("Break : Nan")
			break
		end
	end
end

@views function compute_meca(qDx, qDy,dϕdt, dPedt, Pe, ϕ, k_r, η_f, η_r, Bf, npow, dx, dy, dt, Δρg)
	qDx   .= (k_r/η_f) .* avx(ϕ[:,2:end-1]).^npow .* (.- diff(Pe[:,2:end-1],dims=1) ./ dx)
    qDy   .= (k_r/η_f) .* avy(ϕ[2:end-1,:]).^npow .* (.- diff(Pe[2:end-1,:],dims=2) ./ dy .- Δρg)
    dPedt[2:end-1,2:end-1] .= (1 ./ (ϕ[2:end-1,2:end-1] .* Bf)) .* ( 1 ./ (1 .- ϕ[2:end-1,2:end-1]) .* dϕdt[2:end-1,2:end-1] .- (diff(qDx,dims=1) ./ dx .+ diff(qDy,dims=2) ./ dy))
    Pe[2:end-1,2:end-1] .+= dt .* dPedt[2:end-1,2:end-1]
    dϕdt  .= .-(1 .- ϕ) .* (ϕ ./ η_r) .* Pe 
    ϕ[2:end-1,2:end-1]  .+= dt .* dϕdt[2:end-1,2:end-1]
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

@views function ploting_ϕ(xc,yc,ϕ,qDx,qDy,dϕdt,Pe,dPedt,x_tick,real_time)
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

porowave_1D_pore(true)
