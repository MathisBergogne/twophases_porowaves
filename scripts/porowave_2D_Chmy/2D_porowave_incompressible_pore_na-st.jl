using Chmy, Chmy.Architectures, Chmy.Grids, Chmy.Fields, Chmy.BoundaryConditions, Chmy.GridOperators, Chmy.KernelLaunch
using KernelAbstractions
using Printf
using CairoMakie

@views avx(A) = 0.5 .* (A[1:end-1,:] .+ A[2:end,:])
@views avy(A) = 0.5 .* (A[:,1:end-1] .+ A[:,2:end])

@views function Test_compute_qD()
	lx   = 2500      # longeur du modele en x
	ly   = 10000     # longeur du modele en y
	k_r  = 1e-13     # permeabilite initial m2 / viscosité fluide Pa s
	η_f	 = 1e2       # viscosité fluide Pa s
	η_r  = 1e15      # matrix bulk viscosity
	ϕ_bg = 0.01      # porosité du fond
	loc_ϕ = 0.85	 # localisation of the intitial porosity on the box (0-top/1-bottom)
	w	 = 500		 # shape aspect of the initial melt 
	npow = 3         # exposant de la porosité de l'équation Carman–Kozeny
	ϕA   = 0.1       # différence de porosité entre le fond et une bulle de porosité
	ρs		= 2700							# matrix volumetric mass (kg.m-3)
	ρf		= 2300							# mvolumetric mass fluid (kg.m-3)
	g		= 9.81							# gravity (m.s-2)
	Bf   = 0.1       # compressibilité du fluide 
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
	nt   = 10      # nombre de pas de temps
	dt = 1e-6 * tsc
	nvis = 1 # nt / 100
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
    fig = Figure(; size=(800, 600))
    axs = (ϕ   = Axis(fig[1, 1][1, 1]; aspect=DataAspect(), ylabel="y", xlabel="x", title="ϕ"),
           qDx = Axis(fig[1, 2][1, 1]; aspect=DataAspect(), ylabel="y", xlabel="x", title="qD.x"),
           qDy = Axis(fig[1, 3][1, 1]; aspect=DataAspect(), ylabel="y", xlabel="x", title="qD.y"))
    plt = (ϕ   = heatmap!(axs.ϕ  , xvertices(grid), ycenters(grid), interior(ϕ) |> Array; colormap=:turbo),
           qDx = heatmap!(axs.qDx, xvertices(grid), ycenters(grid), interior(qD.x) |> Array; colormap=:turbo),
           qDy = heatmap!(axs.qDy, xcenters(grid), yvertices(grid), interior(qD.y) |> Array; colormap=:turbo))
    Colorbar(fig[1, 1][1, 2], plt.ϕ)
    Colorbar(fig[1, 2][1, 2], plt.qDx)
    Colorbar(fig[1, 3][1, 2], plt.qDy)
    display(fig)
	#loop
	for it = 1:nt 
		# compute qD.x and qD.y
		launch(arch, grid, compute_qD! => (qD, ϕ, Pe, k_r, η_fg, Δρg, npow, grid))	
		launch(arch, grid, update_porosity! => (ϕ, Pe, dϕdt, η_rg, dt, grid); bc=batch(grid, ϕ => Neumann(); exchange=ϕ))
		launch(arch, grid, update_pressure! => (qD, ϕ, Pe, dϕdt, dPedt, Bf, dt, grid); bc=batch(grid, Pe => Neumann(); exchange=Pe))
	    println("ϕ(max,min) : (",maximum(ϕ),",",minimum(ϕ),")")
		# visualisation
	    if (it % nvis) == 0
			# update visualisation
	    	KernelAbstractions.synchronize(backend)
	    	plt.ϕ[3] = interior(ϕ) |> Array
	    	plt.qDx[3] = interior(qD.x) |> Array
	    	plt.qDy[3] = interior(qD.y) |> Array
	    	display(fig)
		end
	end
end

@kernel inbounds = true function compute_qD!(qD, ϕ, Pe, k_r,η_fg, Δρg, npow, g::StructuredGrid, O)
    I = @index(Global, NTuple)
    I = I + O
    qD.x[I...] = lerp(((k_r[I...]/η_fg[I...]) .* ϕ[I...].^npow),location(qD.x) , g, I...) * ( ∂x(Pe, g, I...) )
    qD.y[I...] = lerp(((k_r[I...]/η_fg[I...]) .* ϕ[I...].^npow),location(qD.y) , g, I...) * ( ∂y(Pe, g, I...) - Δρg)
end

@kernel inbounds = true function update_porosity!(ϕ, Pe, dϕdt, η_rg, dt, g::StructuredGrid, O)
    I = @index(Global, NTuple)
    I = I + O
	dϕdt[I...] = (1 - ϕ[I...])* (ϕ[I...] / η_rg[I...]) * Pe[I...]
	ϕ[I...] += dt * dϕdt[I...]
end

@kernel inbounds = true function update_pressure!(qD, ϕ, Pe, dϕdt, dPedt, Bf, dt, g::StructuredGrid, O)
    I = @index(Global, NTuple)
    I = I + O
	dPedt[I...] = (1 / (ϕ[I...] .* Bf)) * ( 1 / (1 .- ϕ[I...])) * dϕdt[I...] - (∂x(qD.x, g, I...) + (∂y(qD.y, g, I...)))
	Pe[I...] += dt * dPedt[I...]
end

Test_compute_qD()
