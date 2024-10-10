using Plots, Plots.Measures, Printf
default(xmirror = true,size=(900, 600), framestyle=:box, label=false, grid=false, margin=10mm, lw=3, labelfontsize=20, tickfontsize=16, titlefontsize=20)

@views avx(A) = 0.5 .* (A[1:end-1,:] .+ A[2:end,:])
@views avy(A) = 0.5 .* (A[:,1:end-1] .+ A[:,2:end])

@views function porowave_1D_pore(do_visu)
	lx   = 2500       # longeur du modele en x
	ly   = 10000      # longeur du modele en y
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
	Bf   = 0.1       # compressibilité du fluide 
	# independent physics
	lc		= sqrt(k_r / η_f * η_r / ϕ_bg)	# m
	Δρg		= (ρs - ρf) * g					# Pa.m-1
	η_ϕbg	= η_r / ϕ_bg
	psc		= Δρg * lc				# Pa
	tsc		= η_ϕbg / psc	# s
	println("lc : ",lc,";	tsc : ",tsc)
	# numerics
	nx   = 50       # nombre de noeuds de modèles en x
	ny   = 200       # nombre de noeuds de modèles en y
	dx   = lx / nx   # pas d'espace en x
	dy   = ly / ny   # pas d'espace en y
	xc   = LinRange((-lx+dx)/2,(lx-dx)/2,nx)   # coordonnées en x des noeuds du modèles
	yc   = LinRange(-ly+(dy)/2,-dy/2,ny)   # coordonnées en y des noeuds du modèles
	nt   = 1e8      # nombre de pas de temps
	dt = 1e-4 * tsc 
	nvis = nt / 100
	year2sec = 365*24*3600
	### initialisation ###
	# mecanics
	ϕ		= @. ϕ_bg + ϕA * exp(-((xc * loc_ϕ) / w) ^ 2 -((yc' + ly * loc_ϕ) / w) ^ 2)
	Pe     = zeros(nx,ny)
	qDx    = zeros(nx - 1,ny - 2)
	qDy    = zeros(nx - 2,ny - 1)
	dPedt  = zeros(nx - 2,ny - 2)
	dϕdt   = zeros(nx - 2,ny - 2)
	ploting_ϕ(xc,yc,ϕ,Pe,0)
	#loop
	for it = 1:nt 
	    compute_meca(qDx, qDy,dϕdt, dPedt, Pe, ϕ, k_r, η_f, η_r, Bf, npow, dx, dy, dt, Δρg)
	    # visualisation
	    if (it % nvis) == 0 && do_visu || isnan(ϕ[2])
			  ploting_ϕ(xc,yc,ϕ,Pe,dt*it/year2sec) 
	    end
	end

end

@views function compute_meca(qDx, qDy,dϕdt, dPedt, Pe, ϕ, k_r, η_f, η_r, Bf, npow, dx, dy, dt, Δρg)
	qDx   .= (k_r/η_f) .* avx(ϕ[:,2:end-1]).^npow .* (.- diff(Pe[:,2:end-1],dims=1) ./ dx)
    qDy   .= (k_r/η_f) .* avy(ϕ[2:end-1,:]).^npow .* (.- diff(Pe[2:end-1,:],dims=2) ./ dy .- Δρg)
    dϕdt  .= .-(1 .- ϕ[2:end-1,2:end-1]) .* (ϕ[2:end-1,2:end-1] ./ η_r) .* Pe[2:end-1,2:end-1] 
    dPedt .= (1 ./ (ϕ[2:end-1,2:end-1] .* Bf)) .* ( 1 ./ (1 .- ϕ[2:end-1,2:end-1]) .* dϕdt .- (diff(qDx,dims=1) ./ dx .+ diff(qDy,dims=2) ./ dy))
    Pe[2:end-1,2:end-1] .+= dt .* dPedt
    ϕ[2:end-1,2:end-1]  .+= dt .* dϕdt
end

@views function ploting_ϕ(xc,yc,ϕ,Pe,time)
	display(heatmap(xc,yc,ϕ',title="Melt";xlims=(xc[1],xc[end]),ylims=(yc[1],yc[end]),aspect_ratio=1,c=:turbo)) ## ,clim=(0,0.6)
end

porowave_1D_pore(true)
