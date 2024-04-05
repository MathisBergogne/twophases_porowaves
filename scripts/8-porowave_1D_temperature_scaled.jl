using Plots, Plots.Measures, Printf, ElasticArrays

@views avx(A) = 0.5 .* (A[1:end-1] .+ A[2:end])

@views function porowave_1D_speudotime_auto(do_visu)
	#Geology
	#valeurs dimensionnées
	#lx_g	= 5000		# dimmension du modèle (m)
	k_g		= 1e-12		# perméabilité de la matrice (m2)
	ηf_g	= 1e4		# viscosité du fluide (Pa s) 
	ηs_g	= 1e18		# viscosité du solide (Pa s)
	ρs		= 2700		# masse volumique melt (kg m-3)
	ρf		= 2300		# masse volumique rocks (kg m-3)
	g		= -9.81		# gravite terrestre (m s-2)
	λ_g		= 1e-6		# Difusivité thermique (m2 s-1)
	β_g		= 1e-10		# compressibilité du fluide (Pa-1)
	# valeurs sans dimensionsΔρg
	lx_lc	= 100.0	  # 
	Δϕ		= 0.1		# différence de porosité entre le fond et une bulle de porosité
	ϕ_bg	= 1e-2		# porosité du fond
	npow	= 3			# puissnace n pour ϕ

	# independent physics
	lc		= sqrt(k_g / ηf_g * ηs_g / ϕ_bg)	# m
	η_ϕbg	= ηs_g / ϕ_bg						# Pa/m
	Δρg		= (ρs - ρf) * g						# Pa*s
	# scales
	psc		= Δρg * lc
	tsc		= η_ϕbg / psc
	k_ηf0	= lc^2 / η_ϕbg
	#Numerics
	lx		= lx_lc * lc	#lx_g / lc
	λ		= λ_g * lc^2	
	β		= β_g / psc
	kμ0		= k_g / ηf_g		# permeabilite initial m2, viscosité fluide Pa s
	x_pm	= 4/5*lx			# profondeur de la bulle dans la boite	
	Grad_T	= 200 / (4/5*lx)	# gradient géothermique 
	T0		= 750	  
	σ_bulle = lx/50

	# numerics
	nx		= 1000				# nombre de noeuds de modèles en x
	dx		= lx / nx			# pas d'espace
	xc		= LinRange(-dx/2, lx+dx/2,nx)   # coordonnées en x des noeuds du modèles
	nt		= 5e1				# nombre de pas de temps
	dt		= 1e-3 * tsc		# pas de temps
	dτ		= 0.01
	nvis	= 1#nt / 100
	maxiter = 20nx
	ncheck	= ceil(Int,0.05nx)
	ϵtol	= [1e-8, 1e-8]
	re_D	= 4π
	cfl		= 1 / 1.1
	vsdτ	= cfl * dx
	# initialisation
	ϕ		= @. Δϕ  * exp(-(xc - (lx*9/10 + dx / 2))^2 /σ_bulle^2) + ϕ_bg
	φ_init	= copy(ϕ)	
	T		= @. T0 + (Float64(xc) - x_pm) * Grad_T
	T		= @. min(T,T0)
	T_init	= copy(T)
	T_old	= zeros(nx)
	r_T		= zeros(nx - 2)
	dTdt	= zeros(nx - 2)
	qT		= zeros(nx - 1)
	Pe		= zeros(nx)
	Pe_old	= zeros(nx)
	qD		= zeros(nx - 1)
	RPe		= zeros(nx - 2)
	RqD		= zeros(nx - 1)
	η_ϕ		= zeros(nx)
	k_ηf	= zeros(nx)
	lc_loc	= zeros(nx)
	re		= zeros(nx)
	ρ_dτ	= zeros(nx)
	βf_dτ	= zeros(nx)
	η_ϕτ	= zeros(nx)
	k_ηfτ	= zeros(nx - 1)
	re_T	= zeros(nx)
	θ_dτ_T	= zeros(nx)
	β_dτ_T = zeros(nx)
	Count_iter = Float64[]

	#loop
	for it = 1:nt 
		T_old  .= T
		Pe_old .= Pe
		iter = 1; errs = 2 .* ϵtol;iters_evo = Float64[];errs_evo = ElasticArray{Float64}(undef, 2, 0)
		while any(errs .> ϵtol) && iter <= maxiter
			η_ϕ	.= η_ϕbg .* (ϕ_bg ./ ϕ)
			k_ηf   .= k_ηf0 .* (ϕ ./ ϕ_bg) .^ npow

			re_T	= π .+ sqrt.(π^2 .+ lx^2 ./ λ ./ dt)
			θ_dτ_T  = lx ./ re_T ./ cfl ./ dx
			β_dτ_T  = (re_T .* λ) ./ (cfl .*dx .* lx)

			lc_loc .= sqrt.(k_ηf .* η_ϕ)
			re	 .= π .+ sqrt.(π .^ 2 .+ (lx ./ lc_loc) .^ 2)
			ρ_dτ   .= re .* η_ϕ ./ lx ./ dx ./ cfl #ρ_dτ   .= re .* 1.0 ./ (1.0 ./ η_ϕ .+ β ./ dt) ./ lx ./ dx ./ cfl
			βf_dτ  .= 1 ./ vsdτ^2 ./ ρ_dτ
			η_ϕτ   .= 1.0 ./ (βf_dτ .+ 1.0 ./ η_ϕ ) #η_ϕτ   .= 1.0 ./ (βf_dτ .+ 1.0 ./ η_ϕ .+ β ./ dt)
			k_ηfτ  .= 1.0 ./ (avx(ρ_dτ) .+ 1.0 ./ avx(k_ηf))
			
			# fluid pressure update
			RPe   .= .+ diff(qD) ./ dx .- Pe[2:end-1] ./ η_ϕ[2:end-1] 
			Pe[2:end-1] .+= η_ϕτ[2:end-1] .* RPe
			RqD	 .= .-qD ./ avx(k_ηf) .+ diff(Pe) ./ dx .- Δρg
			qD	 .+= k_ηfτ .* RqD
			
			# temperature update
			qT  .-= (qT .+ λ .*(diff(T)./dx))./(1.0 .+ θ_dτ_T)
			dTdt[2:end-1] .= (T[3:end-2] .- T_old[3:end-2])./dt .+ 
				(max.(qD[2:end-2],0.0).*diff(T[2:end-2])./dx .+
				 min.(qD[3:end-1],0.0).*diff(T[3:end-1])./dx )
			T[2:end-1] .-= (dTdt .+ diff(qT)./dx)./(1.0/dt .+ β_dτ_T)
			#T[[1,end]]  .= T[[2,end-1]]
		
			if iter % ncheck == 0
				errs[1] = maximum(abs.(RPe))
				errs[2] = maximum(abs.(RqD))
				append!(errs_evo, errs)
				push!(iters_evo, iter / nx)
				@printf("  iter = %d, iter/nx = %1.3e, errs = [ %1.3e %1.3e ]\n", iter, iter / nx, errs...)
			end
			iter += 1
		end
		push!(Count_iter, iter / nx)
		#porosity update
		ϕ .-= dt * Pe ./ η_ϕ # ϕ[2:end-1]  .+= (1 .- ϕ[2:end-1]) .* dt .* Pe[2:end-1] ./ η_ϕ[2:end-1] .+ ϕ[2:end-1] .* β .* (Pe[2:end-1] .- Pe_old[2:end-1]) ./ dt 

		#condition de break
		if isnan(ϕ[2]) 
			println("break (φ=NaN)")
			break 
		end
		# visualisation 
		if (it % nvis) == 0 && do_visu 
			p1 = plot([ϕ,φ_init], xc, title="Melt",yaxis= :flip, xlims=(0,0.11))
			p2 = plot(Pe, xc, title="Effective Pressure",yaxis= :flip) 
			p3 = plot([T,T_init], xc, title="Temperature",yaxis= :flip, xlims=(0,800))
#			p4 = plot(iters_evo,errs_evo,title="iter / nx",yaxis=:log10, minorgrid=true, marker=:circle)
#			p5 = plot(1:it,Count_iter,title="iter / nx", marker=:circle, ylims=(minimum(Count_iter)*0.9,maximum(Count_iter)*1.1))
			display(plot(p1,p2,p3))
		end
	end
end

porowave_1D_speudotime_auto(true)
