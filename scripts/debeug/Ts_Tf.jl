using Plots, Plots.Measures, Printf, ElasticArrays
default(size=(900, 600), framestyle=:box , label=false,  margin=5mm , lw=2, labelfontsize=20, tickfontsize=10, titlefontsize=16)

@views avx(A) = 0.5 .* (A[1:end-1] .+ A[2:end])

# @views function bubble_heigh(ϕ, lim_bubb)
# 	if isnan(ϕ[2]) 
# 		ϕ_heigh = Missing()
# 	elseif minimum(ϕ) < 0
# 		ϕ_heigh = Missing()
# 	else
# 		ϕ_heigh =  maximum(findall(x->x==maximum(ϕ[1:lim_bubb]), ϕ[1:lim_bubb]))
# 		if ϕ[ϕ_heigh] == ϕ[lim_bubb]
# 			ϕ_heigh =  maximum(findall(x->x==maximum(ϕ), ϕ))
# 		end
# 	end
# 	return ϕ_heigh
# end

@views function bubble_heigh(ϕ,lim_bubb)
	dϕ = ϕ[1:end-1]-ϕ[2:end]
	if isnan(ϕ[2]) 
		ϕ_heigh = Missing()
	elseif minimum(ϕ) < 0
		ϕ_heigh = Missing()
	else
		for k = lim_bubb:-1:2
			if dϕ[k] > 0 && dϕ[k-1] < 0
				ϕ_heigh = k
				break
			end
		end
	end
	return ϕ_heigh
end

@views function break_conditions(ϕ,ϕ_init,xc,qD,T,T_init,real_time,it,nx,lim_bubb,dx,box_top,year2sec)
	ϕ_heigh = bubble_heigh(ϕ, lim_bubb)
	if isnan(ϕ[2]) 
		println("break (φ = NaN at : ",real_time," y; iter : ",it,")")
		return true
	elseif minimum(ϕ) < 0
		neg_ϕ = findall(x->x<0, ϕ) 
		println("break (négative porosity located at ", - neg_ϕ .* dx .- box_top , "; at : ",real_time," y; iter : ",it,")")
		return true 
	elseif ϕ_heigh >= lim_bubb 
		println("break (bublle rise the top at : ",real_time," y; iter : ",it,")")
		println(@sprintf("ϕ_heigh : %f; nx : %f",ϕ_heigh,nx))
		ΔE = (sum(T)-sum(T_init))/sum(T_init)*100
		Δϕ = (sum(ϕ)-sum(ϕ_init))/sum(ϕ_init)*100
		dϕ = ϕ[1:end-1]-ϕ[2:end]
		title_plot = @sprintf(" time : %1.1e y; ΔE : %1.1e %s; Δϕ : %1.1e %s",real_time/year2sec,ΔE,"%",Δϕ,"%")
		p1 = plot(ϕ.*100, xc, title="% of Melt", xlims=(0,30))
		p1 = scatter!([ϕ[ϕ_heigh].*100], [xc[ϕ_heigh]])
		p5 = plot(qD, avx( xc), title="qD")
		p6 = plot(dϕ, avx( xc), title="dϕ")
		display(plot(p1,p5,p6,plot_title=title_plot; layout=(1, 3)))
		return true
	#  elseif T[ϕ_heigh] <= T_init[ϕ_heigh] && (nx - ϕ_heigh) <= lim_bubb 
	# 	 println("break (bublle cooldown at : ",real_time," )")
	# 	 return true
	end
	return false
end

@views function ploting(ϕ,ϕ_init,Ts,Tf,Ts_init,Tf_init,qD,ηf_r,η_r,xc,real_time,ϕ_heigh,lim_bubb,dx,lx,box_top,it,year2sec,cte_η_f,cte_η_r,save_plots)
	ΔE = (sum(Ts)-sum(Ts_init))/sum(Ts_init)*100 + (sum(Tf)-sum(Tf_init))/sum(Tf_init)*100
	Δϕ = (sum(ϕ)-sum(ϕ_init))/sum(ϕ_init)*100
	dϕ = ϕ[1:end-1]-ϕ[2:end]
	title_plot = @sprintf(" time : %1.1e y; ΔE : %1.1e %s; Δϕ : %1.1e %s",real_time/year2sec,ΔE,"%",Δϕ,"%")
	# set up plots
	p1 = plot([ϕ,ϕ_init].*100, xc, title="% of Melt", xlims=(0,30))
	if do_p_bb_h p1 = scatter!([ϕ[ϕ_heigh].*100], [xc[ϕ_heigh]]) end
	if cte_η_f == false p2 = plot(ηf_r, xc, title="Fluid viscosity", xaxis=:log10) end
	if cte_η_r == false p3 = plot(η_r, xc, title="Matrix viscosity", xaxis=:log10) end
	p4 = plot([Ts Ts_init], xc, title="Matrix Temperature")
	if do_p_bb_h p4 = scatter!([Ts[ϕ_heigh]], [xc[ϕ_heigh]]) end
	p5 = plot([Tf Tf_init], xc, title="Fluid Temperature")
	if do_p_bb_h p4 = scatter!([Tf[ϕ_heigh]], [xc[ϕ_heigh]]) end
	p6 = plot(qD, avx(xc), title="qD")
	# Display
	if cte_η_f && cte_η_r 
		#display(plot(p1,p4,plot_title=title_plot; layout=(1, 2)))
		display(plot(p1,p4,p5,plot_title=title_plot; layout=(1, 3)))
		if save_plots savefig(@sprintf("Constant_viscous_%04d",it)) end
	elseif cte_η_f 
		display(plot(p1,p3,p4,plot_title=title_plot; layout=(1, 3)))
		if save_plots savefig(@sprintf("Frank_rock_viscous_%04d",it)) end
	elseif cte_η_r 
		display(plot(p1,p2,p4,plot_title=title_plot; layout=(1, 3)))
		if save_plots savefig(@sprintf("Frank_Fluid_viscous_%04d",it)) end
	else 
		display(plot(p1,p2,p3,p4,plot_title=title_plot; layout=(1, 4)))
		#display(plot(p1,p4,p6,plot_title=title_plot; layout=(1, 3)))
		if save_plots savefig(@sprintf("Frank_viscous_%04d",it)) end
	end
end

@views function Temperature_E(do_visu, do_p_itr, do_p_bb_h, do_p_vel,save_plots,cte_η_f,cte_η_r, perm, η0_r, θ_r, η0_f, θ_f)
	# real values
	lx		= 5e3							# lenght of the model (m)
	box_top = 5e3
	ρs		= 2700							# matrix volumetric mass (kg.m-3)
	ρf		= 2300							# fluid volumetric mass (kg.m-3)
	Cps		= 1000							# matrix Heat capacity (J K-1 kg-1).
	Cpf		= 4000 	### Cps*ρs/ρf			# fluid Heat capacity (J K-1 kg-1).
	ρCpf 	= ρf .* Cpf
	ρCps 	= ρs .* Cps
	λs		= 3								# matrix Thermal conductivity (kg m s-3 K-1) 
	λf		= 2								# fluid Thermal conductivity (kg m⋅s-3 K-1)
	β_f		= 1e-10							# compressibilité du fluide (Pa-1)
	g		= 9.81							# gravity (m.s-2)
	w		= 100							# shape aspect of the initial melt 
    ϕ_bg	= 1e-2							# back groud amount of melt 
	ϕA		= 0.2							# amount of initial melt
	npow	= 3								# ...
	year2sec = 365*24*3600
	loc_ϕ	= 0.85
	# numerics thermal gradiant
	loc_T	= 4/5
	T_top	= 500
	T0		= 750
	x_pm	= loc_T*lx			# profondeur de la bulle dans la boite	
	Grad_T	= 0 # (T0 - T_top) / (loc_T*lx)	# gradient géothermique 	
	# melt and rock properties
	k_r		= perm							# matrix permeability (m2)
	ηf_min	= η0_f .* exp.(- θ_f .* T0)	# fluid viscosity (Pa.s) 
	ηr_min	= η0_r .* exp.(- θ_r .* T0)	# matrix viscosity (Pa.s)
	# independent physics
	lc		= sqrt(k_r / ηf_min * ηr_min / ϕ_bg)	# m
	Δρg		= (ρs - ρf) * g					# Pa.m-1
	η_ϕbg_min	= ηr_min / ϕ_bg
	# numerics space
	dx_approx = lc
	nx		= max(Int(round(lx/dx_approx))+1,200)
	dx		= lx/(nx-1)
	xc		= LinRange( - (lx + box_top + dx / 2), - box_top + dx / 2, nx)
	# init Temperature
	Ts 		= Tf		= zeros(nx)
	Ts 		= Tf		= @. T0 + (Float64(-xc - box_top) - x_pm) * Grad_T
	Ts ,Tf	= @. min(Ts,T0), @. min(Ts,T0)
	Ts_init	= copy(Ts)
	Tf_init	= copy(Tf)
	Ts_old 	= Tf_old	= zeros(nx)
	qd_s 	= qd_f		= zeros(nx - 1)
	qa_s 	= qa_f		= zeros(nx - 1)
	# scales
	Lsc     = 1.					# m
	psc		= Δρg * lc				# Pa
	tsc		= η_ϕbg_min / psc	# s
	# numerics time
	nt 		= 1
	dt_meca = 1e-3 * tsc
	nvis	= 1  
	real_time = 0.0
	dt_diff	= 0.0
	dt_adv	= 0.0
	dt		= 0.0
	if do_p_vel dt_vel	= 0.0 end
	# numerics other
	maxiter	= 50nx
	ϵtol	= [1e-8, 1e-8]
	cfl		= 1 / 1.1
	ncheck	= ceil(Int, 0.25nx)
	# preprocessing
	vsdτ	= cfl * dx
	# init porosity
	ϕ		= @. ϕ_bg + ϕA * exp(-((xc + box_top + lx * loc_ϕ) / w) ^ 2)
	#ϕ[1] 	= 1e-7
	ϕ_init	= copy(ϕ)
	ϕ_vel 	= Float64[] # ElasticArray{Float64}(undef, 1, 0)
	qD_vel	= Float64[]
	lim_bubb = Int(round(0.75 * nx))
	ϕ_heigh = bubble_heigh(ϕ, lim_bubb)
	# init pressure
	Pe		= zeros(nx)
	qD		= zeros(nx - 1)
	# init speudo transient
	η_ϕ		= zeros(nx)
	k_ηf	= zeros(nx)
	RPe		= zeros(nx - 2)
	RqD		= zeros(nx - 1)
	η_ϕτ	= zeros(nx)
	βf_dτ	= zeros(nx)
	k_ηfτ	= zeros(nx - 1)
	ρ_dτ	= zeros(nx)
	re		= zeros(nx)
	lc_loc	= zeros(nx)
	save_dt = zeros(nt)
	# action
	for it = 1:nt
		iter = 1
		errs = 2.0 .* ϵtol
		errs_evo = ElasticArray{Float64}(undef, 2, 0)
		iters_evo = Float64[]
		# Ts_old, Tf_old  .= Ts, Tf
		# prorosity dependance
		ρCpT = ρf .* Cpf .* ϕ .* (1 .- ϕ) .+ ρs .* Cps
		λT = λf .* avx(ϕ) .+ λs .* (1 .- avx(ϕ))
		χT = λT ./ avx(ρCpT)
		#test eta Frank
		η_r		= η0_r .* exp.(- θ_r .* Ts)	# matrix viscosity (Pa.s)
	    ηf_r	= η0_f .* exp.(- θ_f .* Tf)	# fluid viscosity (Pa.s)	
		η_ϕbg	= η_r / ϕ_bg					# Pa.s	
		k_ηf0	= k_r ./ ηf_r					# m3.s2.kg-1	
		# time step 
		dt_diff		= 0.5*dx^2 / maximum(χT) 
		dt_adv		= 0.1*dx/maximum(abs.(qD))
		dt = min(dt_adv,dt_meca) # dt_diff #
		real_time   += dt 
		save_dt[it] = dt
		if do_p_vel dt_vel += dt end
		# speudo transient solve
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
			if melt_input == false; qD[1] 	 = 0.0 end
			# check error
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
			# display(plot(Ts, xc, title="Ts"),plot(Tf, xc, title="Tf"))
			qd_s .= .- λs .* (diff(Ts) ./ (dx)) 
			qd_f .= .- λf .* (diff(Tf) ./ (dx))
			qa_s .= qD .* avx(Ts) .* ρCps
			qa_f .= .- qD .* avx(Tf) .* ρCpf
			Ts[2:end-1] .= Ts[2:end-1] - diff(qa_s)./dx .* dt_diff ./ (ρCps .* ϕ[2:end-1])
			Tf[2:end-1] .= Tf[2:end-1] - diff(qa_f)./dx .* dt_diff ./ (ρCpf .* (1 .- ϕ[2:end-1]))  
			# Ts[1:Int(round((1-loc_T)*nx))] .= Tf[1:Int(round((1-loc_T)*nx))] .= T0
		end
		#porosity update
		ϕ .-= dt * (1 .- ϕ) .* Pe ./ η_ϕ 
		if melt_input == false; ϕ[1]	= 1e-4 end
		#Break ?
		# do_break = break_conditions(ϕ,ϕ_init,xc,qD,T,T_init,real_time/year2sec,it,nx,lim_bubb,dx,box_top,year2sec)
		# if do_break
			# nt = it
			# println("viscosity of the matrix : ",η_r[ϕ_heigh],";	viscosity of the melt : ", ηf_r[ϕ_heigh])
			# ploting(ϕ,ϕ_init,Ts,Tf,Ts_init,Tf_init,qD,ηf_r,η_r,xc,real_time,ϕ_heigh,lim_bubb,dx,lx,box_top,it,year2sec,cte_η_f,cte_η_r,save_plots)
			# break
		# end
		# visualisation 
		if (it % nvis) == 0 && do_visu
			old_ϕ_heigh = ϕ_heigh
			ϕ_heigh = bubble_heigh(ϕ, lim_bubb)
			if do_p_vel 
				push!(ϕ_vel, (ϕ_heigh - old_ϕ_heigh) * dx / (dt_vel / year2sec)) 
				push!(qD_vel, qD[ϕ_heigh])
				dt_vel = 0.0
			end		
			ploting(ϕ,ϕ_init,Ts,Tf,Ts_init,Tf_init,qD,ηf_r,η_r,xc,real_time,ϕ_heigh,lim_bubb,dx,lx,box_top,it,year2sec,cte_η_f,cte_η_r,save_plots)
		end
	end
	if do_p_vel 
		display(plot(ϕ_vel, nvis:nvis:nt, xlims=(min(-0.005,minimum(ϕ_vel)*-0.9),maximum(ϕ_vel)*1.1)))
	end
end

# some boolean
do_visu		= true		# to make plot
do_p_bb_h	= true		# to visualize the bubble on plot
save_plots	= false		# to save plots
do_p_vel	= false		# to plot the velocity over time
do_p_itr	= false 	# to print iterations
melt_input	= true
#set constant viscosity ?
cte_η_f = true			# to set a contant fluid viscosity
cte_η_r = true			# to set a contant rock viscosity
# set rock viscosity parameters
if cte_η_r
	η0_r 	= 1e17 		# constant rock viscosity 
	θ_r 	= 0.0		# do not change
else 
	η0_r 	= 1e30		# viscosity at 0°C
	θ_r 	= 0.038		# coefficient with temperature in exponential 
end 
# set fluid viscosity parameters
if cte_η_f
	η0_f	= 1e5 		# constant fluid viscosity
	θ_f 	= 0.0		# do not change
else 
	η0_f	= 1e13		# viscosity at 0°C
	θ_f 	= 0.03		# coefficient with temperature in exponential 
end 

perm 	= 1e-11 		# matrix permeabiblity

Temperature_E(do_visu, do_p_itr, do_p_bb_h, do_p_vel,save_plots,cte_η_f,cte_η_r, perm, η0_r, θ_r, η0_f, θ_f)
