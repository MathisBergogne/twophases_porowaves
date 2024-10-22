using Chmy, Chmy.Architectures, Chmy.Grids, Chmy.Fields, Chmy.BoundaryConditions, Chmy.GridOperators, Chmy.KernelLaunch
using KernelAbstractions
using Printf
using CairoMakie

#using Plots, Plots.Measures, Printf
#default(xmirror = true,size=(900, 600), framestyle=:box, label=false, grid=false, margin=10mm, lw=3, labelfontsize=20, tickfontsize=16, titlefontsize=20)

@views avx(A) = 0.5 .* (A[1:end-1,:] .+ A[2:end,:])
@views avy(A) = 0.5 .* (A[:,1:end-1] .+ A[:,2:end])

@views function Test_compute_qD_Chmy()
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
	dt = 1e-4 * tsc
	### initialisation ###
	# mecanics
	ϕ		= @. ϕ_bg + ϕA * exp(-((xc * loc_ϕ) / w) ^ 2 -((yc' + ly * loc_ϕ) / w) ^ 2)
	Pe     = zeros(nx,ny)
	qDx    = zeros(nx - 1,ny - 2)
	qDy    = zeros(nx - 2,ny - 1)
	# Chmy param
	backend=CPU()
	arch = Arch(backend)
	nxy=(nx, ny)
    # geometry
    grid   = UniformGrid(arch; origin=(-lx/2, -ly), extent=(lx, ly), dims=nxy)
    launch = Launcher(arch, grid; outer_width=(16, 8))
	# fields
	ϕ_Chmy = Field(backend, grid, Center())
	Pe_Chmy = Field(backend, grid, Center())
	k_r_Chmy = Field(backend, grid, Center())
	η_f_Chmy = Field(backend, grid, Center())
    qD = VectorField(backend, grid)
    # initial conditions
	set!(ϕ_Chmy, grid, (x, y) -> ϕ_bg + ϕA * exp(-((x * loc_ϕ) / w) ^ 2 -((y + ly * loc_ϕ) / w) ^ 2))
	bc!(arch, grid, ϕ_Chmy => Neumann(); exchange=ϕ_Chmy)
    set!(Pe_Chmy, grid, (_, _) -> rand())
	bc!(arch, grid, Pe_Chmy => Neumann(); exchange=Pe_Chmy)
    set!(k_r_Chmy, grid, (_, _) -> k_r)
	bc!(arch, grid, k_r_Chmy => Neumann(); exchange=k_r_Chmy)
    set!(η_f_Chmy, grid, (_, _) -> η_f)
	bc!(arch, grid, η_f_Chmy => Neumann(); exchange=η_f_Chmy)
	# visualisation
    fig = Figure(; size=(800, 600))
    axs = (ϕ   = Axis(fig[1, 1][1, 1]; aspect=DataAspect(), ylabel="y", xlabel="x", title="ϕ"),
           ϕ_C = Axis(fig[2, 1][1, 1]; aspect=DataAspect(), ylabel="y", xlabel="x", title="ϕ_C"),
           qDx = Axis(fig[1, 2][1, 1]; aspect=DataAspect(), ylabel="y", xlabel="x", title="qDx"),
           qDX = Axis(fig[2, 2][1, 1]; aspect=DataAspect(), ylabel="y", xlabel="x", title="qD.x"),
           qDy = Axis(fig[1, 3][1, 1]; aspect=DataAspect(), ylabel="y", xlabel="x", title="qDy"),
           qDY = Axis(fig[2, 3][1, 1]; aspect=DataAspect(), ylabel="y", xlabel="x", title="qD.y"))
    plt = (ϕ   = heatmap!(axs.ϕ  ,xc,yc,ϕ;colormap=:turbo),
           ϕ_C = heatmap!(axs.ϕ_C, xvertices(grid), ycenters(grid), interior(ϕ_Chmy) |> Array; colormap=:turbo),
           qDx = heatmap!(axs.qDx,xc[1:end-1],yc[2:end-1],qDx;colormap=:turbo),
           qDX = heatmap!(axs.qDX, xvertices(grid), ycenters(grid), interior(qD.x) |> Array; colormap=:turbo),
           qDy = heatmap!(axs.qDy,xc[2:end-1],yc[1:end-1],qDy;colormap=:turbo),
           qDY = heatmap!(axs.qDY, xcenters(grid), yvertices(grid), interior(qD.y) |> Array; colormap=:turbo))
    Colorbar(fig[1, 1][1, 2], plt.ϕ)
    Colorbar(fig[2, 1][1, 2], plt.ϕ_C)
    Colorbar(fig[1, 2][1, 2], plt.qDx)
    Colorbar(fig[2, 2][1, 2], plt.qDX)
    Colorbar(fig[1, 3][1, 2], plt.qDy)
    Colorbar(fig[2, 3][1, 2], plt.qDY)
    display(fig)
	# compute qDx and qDy
	compute_qD!(qDx, qDy, Pe, ϕ, k_r, η_f, npow, dx, dy, Δρg)
	# compute qD.x and qD.y
	launch(arch, grid, compute_qD_Chmy! => (qD, ϕ_Chmy, Pe_Chmy, k_r_Chmy, η_f_Chmy, Δρg, npow, grid))
	# update visualisation
    KernelAbstractions.synchronize(backend)
    plt.ϕ[3]   = ϕ |> Array
    plt.ϕ_C[3] = interior(ϕ_Chmy) |> Array
    plt.qDx[3] = qDx |> Array
    plt.qDX[3] = interior(qD.x) |> Array
    plt.qDy[3] = qDy |> Array
    plt.qDY[3] = interior(qD.y) |> Array
    display(fig)
end

@kernel inbounds = true function compute_qD_Chmy!(qD, ϕ_Chmy, Pe_Chmy, k_r_Chmy,η_f_Chmy, Δρg, npow, g::StructuredGrid, O)
    I = @index(Global, NTuple)
    I = I + O
    qD.x[I...] = (k_r_Chmy[I...]/η_f_Chmy[I...]) .* ϕ_Chmy[I...].^npow * ( ∂x(Pe_Chmy, g, I...) )
    qD.y[I...] = (k_r_Chmy[I...]/η_f_Chmy[I...]) .* ϕ_Chmy[I...].^npow * ( ∂y(Pe_Chmy, g, I...) - Δρg)
end


@views function compute_qD!(qDx, qDy, Pe, ϕ, k_r, η_f, npow, dx, dy, Δρg)
	qDx   .= (k_r/η_f) .* avx(ϕ[:,2:end-1]).^npow .* (.- diff(Pe[:,2:end-1],dims=1) ./ dx)
    qDy   .= (k_r/η_f) .* avy(ϕ[2:end-1,:]).^npow .* (.- diff(Pe[2:end-1,:],dims=2) ./ dy .- Δρg)
end


Test_compute_qD_Chmy()
