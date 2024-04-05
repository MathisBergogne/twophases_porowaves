using Plots, Plots.Measures, Printf, ElasticArrays
default(size=(900, 600), framestyle=:box , label=false,  margin=5mm , lw=2, labelfontsize=20, tickfontsize=10, titlefontsize=16)

@views avx(A) = 0.5 .* (A[1:end-1] .+ A[2:end])
	
@views function Temperature_E(do_visu, do_p_itr, perm, visco_f, visco_r)
	# input
	k_r		= perm							# matrix permeability (m2)
	ηf_r	= visco_f						# fluid viscosity (Pa.s) 
	η_r		= visco_r						# matrix viscosity (Pa.s)
	# real values
	lx		= 5e3 # 1e4						# lenght of the model (m)
	box_top = 5e3
	ρs		= 2700							# matrix volumetric mass (kg.m-3)
	ρf		= 2300							# mvolumetric mass fluid (kg.m-3)
	Cps		= 1000
	Cpf		= 800
	λs		= 3
	λf		= 2
	g		= 9.81							# gravity (m.s-2)
	w		= 100							# shape aspect of the initial melt 
    ϕ_bg	= 1e-2							# back groud amount of melt 
	ϕA		= 0.1							# amount of initial melt
	npow	= 3								# ...
	k_ηf0	= k_r / ηf_r					# m3.s2.kg-1
	# independent physics
	lc		= sqrt(k_r / ηf_r * η_r / ϕ_bg)	# m
	Δρg		= (ρs - ρf) * g					# Pa.m-1
    η_ϕbg	= η_r / ϕ_bg					# Pa.s
	# scales
	Lsc     = 1.                             # m
	psc		= Δρg * lc						# Pa
	tsc		= η_ϕbg / psc					# s
	χsc		= Lsc^2/tsc
	# numerics
	nt 		= 50
	dt_meca = 1e-3 * tsc
	nvis	= 1 # nt/50  
	real_time = 0.0
	dt_diff	= 0.0
	dt_adv	= 0.0
	dt		= 0.0
	# numerics thermal gradiant
	T_top	= 500
	T0		= 750	
	x_pm	= 4/5*lx			# profondeur de la bulle dans la boite	
	Grad_T	= (T0 - T_top) / (4/5*lx)	# gradient géothermique 
	# numerics time
	dx_approx = lc
	nx		= max(Int(round(lx/dx_approx))+1,200)
	dx		= lx/(nx-1)
	# numerics other
	maxiter	= 50nx
	ϵtol	= [1e-8, 1e-8]
	cfl		= 1 / 1.1
	ncheck	= ceil(Int, 0.25nx)
	# preprocessing
	xc		= LinRange( - (lx + box_top + dx / 2), - box_top + dx / 2, nx)
	vsdτ	= cfl * dx
	# init porosity
	ϕ		= @. ϕ_bg + ϕA * exp(-((xc + box_top + lx * 0.85) / w) ^ 2)
	ϕ_init	= copy(ϕ)
	η_ϕ		= zeros(nx)
	k_ηf	= zeros(nx)
	# init pressure
	Pe		= zeros(nx)
	qD		= zeros(nx - 1)
	RPe		= zeros(nx - 2)
	RqD		= zeros(nx - 1)
	# init speudo transient
	η_ϕτ	= zeros(nx)
	βf_dτ	= zeros(nx)
	k_ηfτ	= zeros(nx - 1)
	ρ_dτ	= zeros(nx)
	re		= zeros(nx)
	lc_loc	= zeros(nx)
	# init Temperature
	T		= zeros(nx)
	T		= avx(@. T0 + (Float64(-xc - box_top) - x_pm) * Grad_T)
	T	  	= @. min(T,T0)
	T_init 	= copy(T)
	T_old  	= zeros(nx - 1)
	qd_T	= zeros(nx - 2)
	qa_f	= zeros(nx - 2)
	# action
	for it = 1:nt
		iter = 1
		errs = 2.0 .* ϵtol
		errs_evo = ElasticArray{Float64}(undef, 2, 0)
		iters_evo = Float64[]
		T_old  .= T
		# prorosity dependance
		ρCpT = ρf .* Cpf .* avx(ϕ) .+ ρs .* Cps .* (1 .-avx(ϕ))
		ρCpf = (ρf .* Cpf)     
		λT = λf .* ϕ .+ λs .* (1 .- ϕ)
		χT = avx(λT) ./ ρCpT
		# time step 
		dt_diff		= 0.5*dx^2 / maximum(χT) 
		dt_adv		= 0.1*dx/maximum(abs.(qD))
		dt = min(dt_adv,dt_meca) # dt_diff #
		real_time   += dt
		while any(errs .> ϵtol) && iter <= maxiter
			# material properties
			η_ϕ		.= η_ϕbg .* (ϕ_bg ./ ϕ)
			k_ηf	.= k_ηf0 .* (ϕ ./ ϕ_bg) .^ npow
			# numerical parameter update
			lc_loc	.= sqrt.(k_ηf .* η_ϕ)
			re		.= π .+ sqrt.(π .^ 2 .+ (lx ./ lc_loc) .^ 2)
			ρ_dτ	.= re .* η_ϕ ./ lx ./ dx ./ cfl
			βf_dτ	.= 1 / vsdτ^2 ./ ρ_dτ
			η_ϕτ	.= 1.0 ./ (βf_dτ .+ 1.0 ./ η_ϕ)
			k_ηfτ	.= 1.0 ./ (avx(ρ_dτ) .+ 1.0 ./ avx(k_ηf))
			# update of physical fields
			RPe		.= .-Pe[2:end-1] ./ η_ϕ[2:end-1] .- diff(qD) ./ dx
			Pe[2:end-1] .+= RPe .* η_ϕτ[2:end-1]
			RqD		.= .-qD ./ avx(k_ηf) .- diff(Pe) ./ dx .- Δρg
			qD	   .+= RqD .* k_ηfτ
			if iter % ncheck == 0
				errs[1] = maximum(abs.(RPe))
				errs[2] = maximum(abs.(RqD))
				append!(errs_evo, errs)
				push!(iters_evo, iter / nx)
				if do_p_itr @printf("it = %d, iter = %d, iter/nx = %1.3e, errs = [ %1.3e %1.3e ]\n",it, iter, iter / nx, errs...) end
			end
			iter += 1
		end
		# compute themperature
		for i = 1:round(dt/dt_diff)
			# temperature update based on flux
			qd_T .=  .- λT[2:end-1] .* (diff(T) ./ (dx)) 
			qa_f .= ρCpf .* (min.(-qD[1:end-1],0.0) .*T_old[1:end-1].+ max.(-qD[2:end],0.0).*T_old[2:end])
			T[2:end-1] .= T[2:end-1] - diff(qd_T .+ qa_f)./dx .* dt_diff ./ ρCpT[2:end-1] 
			T[1:Int(round((1/5)*nx))] .= 750
		end
		#porosity update
		ϕ .-= dt * Pe ./ η_ϕ
		# visualisation 
		if (it % nvis) == 0 && do_visu 
			p1 = plot([ϕ,ϕ_init].*100, xc, title="% of Melt", xlims=(0,30))
			#p2 = plot(Pe, xc, title="Effective Pressure") 
			p3 = plot([T T_init], avx(xc), title="Temperature", xlims=(450,850))
			display(plot(p1,p3,plot_title=@sprintf(" time : %1.1e y; ΔE : %1.1e %s",real_time/(365*24*3600),sum((T_init.-T)./T_init)/100,"%"); layout=(1, 2)))
			#savefig(@sprintf("Temp_meca_low_viscous%04d",it))
		end
	end
end

do_visu   = true
do_p_itr  = false 

perm 	= 1e-11
visco_f = 1e4
visco_r = 1e17

Temperature_E(do_visu, do_p_itr, perm, visco_f, visco_r)
