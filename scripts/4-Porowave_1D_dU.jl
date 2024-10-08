using Plots, Plots.Measures, Printf, ElasticArrays
default(size=(900, 600), framestyle=:box , label=false,  margin=5mm , lw=2, labelfontsize=20, tickfontsize=10, titlefontsize=16)

###################
#### Functions ####
###################

@views avx(A) = 0.5 .* (A[1:end-1] .+ A[2:end])

@views function results_file(f_result, do_save)
	# !!! If you change " to_save_params " change also " to_save " in function " scaled_porowaves " !!!
	to_save_params =["porosity heigh", "bubble Temperature (ΔT, °C)", "Permeability (m^2)","rock viscosity at 750°C(Pa.s)", "rock viscosity dependance of temperature", "fluid viscosity at 750°C(Pa.s)", "fluid viscosity dependance of temperature",	"dt (s)","nt","dx (m)","nx","lc (m)","Δρg (Pa.m-1)","η_ϕbg (Pa.s)","psc (Pa)","tsc (s)"]
	if isfile(f_result) 
		do_save = false
	else
		touch(f_result)
		CSV.write(f_result, [], writeheader=true, header=to_save_params)
	end 
end

@views function saving_data(to_add, file_name)
	CSV.write(file_name, Tables.table(to_add'), append = true, writeheader = false)
end

# @views function bubble_heigh(ϕ, lim_bubb)
# 	if isnan(ϕ[2]) 
# 		ϕ_heigh = -1
# 	elseif minimum(ϕ) < 0
# 		ϕ_heigh = -1
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
	if isnan(ϕ[2]) || minimum(ϕ) < 0
		ϕ_heigh = -1
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

@views function bubble_heat(T, T_init, ϕ_heigh,ϕ)
	if isnan(ϕ[2]) || minimum(ϕ) < 0
		T_diff = -1
	else
		T_diff = T[ϕ_heigh]	- T_init[ϕ_heigh]
	end
	return T_diff
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

@views function show_caracteristic(Lc, Tsc, xc, dir_plot,save_plots)
	p1 = plot(Lc, xc, title="Caracteristic Lenght")
	exponant = LinRange(log10(round(minimum(Tsc), sigdigits=1, RoundDown)), log10(round(maximum(Tsc), sigdigits=1, RoundUp)), max(2,Int(round(log10(maximum(Tsc)),digits=0))-Int(round(log10(minimum(Tsc)),digits=0))+1))
	xtick_values = []
	xtick_labels = []
	for n in eachindex(exponant)
		push!(xtick_values,round(10^exponant[n])) 
		push!(xtick_labels, @sprintf("%1.0e",xtick_values[n]))
	end
	p2 = plot(Tsc, xc, title="Caracteristic Time", xaxis=:log10, xticks=(xtick_values, xtick_labels), xlims=(round(minimum(Tsc), sigdigits=1, RoundDown),round(maximum(Tsc), sigdigits=1, RoundUp)))
	# Display
		display(plot(p1,p2; layout=(1, 2)))
		if save_plots savefig(@sprintf("%s/Lc_Tsc",dir_plot)) end
end

@views function ploting(ϕ,ϕ_init,T,T_init,ηf_r,η_r,xc,real_time,ϕ_heigh,it,year2sec,save_plots,do_p_bb_h,dir_plot)
	ΔE = (sum(T)-sum(T_init))/sum(T_init)*100
	Δϕ = (sum(ϕ)-sum(ϕ_init))/sum(ϕ_init)*100
	title_plot = @sprintf(" time : %1.1e y; ΔE : %1.1e %s; Δϕ : %1.1e %s",real_time/year2sec,ΔE,"%",Δϕ,"%")
	# set up plots
	p1 = plot([ϕ,ϕ_init].*100, xc, title="% of Melt", xlims=(0,30))
	if do_p_bb_h p1 = scatter!([ϕ[ϕ_heigh].*100], [xc[ϕ_heigh]]) end
	p2 = plot(ηf_r, xc, title="Fluid viscosity", xaxis=:log10)
	p3 = plot(η_r, xc, title="Matrix viscosity", xaxis=:log10)
	p4 = plot([T T_init], xc, title="Temperature")
	if do_p_bb_h p4 = scatter!([T[ϕ_heigh]], [xc[ϕ_heigh]]) end
	# Display
		display(plot(p1,p4,p2,p3,plot_title=title_plot; layout=(1, 4)))
		if save_plots savefig(@sprintf("%s/%04d",dir_plot,it)) end
end

@views function Compute_T(T, dU_dt, qD, ρCpT, λT, ρCpf, Grad_T, loc_T, T0, nx, dx, dt, dt_diff)
	for i = 0:round(dt/dt_diff)
		dU_dt .= ( .- avx(qD) .* ρCpf .* (T[3:end].-T[1:end-2])./(2*dx)
			# ) .- (diff(λT)./dx .* diff(avx(T))./dx
			# ) .- diff( λT .* diff(T)./dx)./dx
			) .- (diff(λT)./dx .* diff(avx(T))./dx
			) .- avx(λT) .* diff(diff(T)./dx)./dx 
		T[2:end-1] .= T[2:end-1] .- ( dU_dt ) .* dt_diff ./ ρCpT[2:end-1] 
		if Grad_T != 0; T[1:Int(round((1-loc_T)*nx))] .= T0; end
	end
end

@views function Compute_meca(ϕ, Pe, dPedt, qD, dϕdt, Δρg,k_ηf, η_r, β_f, npow, dx, dt, melt_input)
	#qD   .= avx(k_ηf .* ϕ.^npow) .* (.- diff(Pe) ./ dx .- Δρg)
    dϕdt  .= .-(1 .- ϕ[2:end-1]) .* (ϕ[2:end-1] ./ η_r[2:end-1]) .* Pe[2:end-1] 
    dPedt .= (1 ./ (ϕ[2:end-1] .* β_f)) .* ( 1 ./ (1 .- ϕ[2:end-1]) .* dϕdt .- diff(qD) ./ dx)
    Pe[2:end-1] .+= dt .* dPedt
    ϕ[2:end-1]  .+= dt .* dϕdt
	if melt_input == false; qD[1] 	 = 0.0; ϕ[1]	= 1e-4 end
end

@views function dU_1D(η0_r, θ_r, η0_f, θ_f, perm, do_visu , do_p_bb_h, save_plots, do_save, do_p_itr, melt_input, dir_plot, f_result)
	# real values
	lx		= 5e3							# lenght of the model (m)
	box_top = 0.0
	ρs		= 2700							# matrix volumetric mass (kg.m-3)
	ρf		= 2300							# mvolumetric mass fluid (kg.m-3)
	Cps		= 790
	Cpf		= 4185 # Cps*ρs/ρf
	λs		= 3
	λf		= 2
	β_f		= 1e-10							# compressibilité du fluide (Pa-1)
	g		= 9.81							# gravity (m.s-2)
	w		= 100							# shape aspect of the initial melt 
    ϕ_bg	= 1e-2							# back groud amount of melt 
	ϕA		= 0.1							# amount of initial melt
	npow	= 3								# ...
	year2sec = 365*24*3600
	loc_ϕ	= 0.85
	# numerics thermal gradiant
	loc_T	= 4/5
	T_top	= 500
	T0		= 750
	x_pm	= loc_T*lx			# profondeur de la bulle dans la boite	
	Grad_T	= (T0 - T_top) / (loc_T*lx)	# gradient géothermique 	
	# melt and rock properties
	k_r		= perm
	# independent physics
	lc		= sqrt(k_r / η0_f * η0_r / ϕ_bg)	# m
	Δρg		= (ρs - ρf) * g					# Pa.m-1
	η_ϕbg_min	= η0_r / ϕ_bg
	# numerics space
	dx_approx = lc
	nx		= max(Int(round(lx/dx_approx))+1,200)
	dx		= lx/(nx-1)
	xc		= LinRange( - (lx + box_top + dx / 2), - box_top + dx / 2, nx)
	# init Temperature
	T		= zeros(nx)
	T		= @. T0 + (Float64(-xc - box_top) - x_pm) * Grad_T
	T	  	= @. min(T,T0)
	T_init 	= copy(T)
	dU_dt	= zeros(nx-2)
	# scales
	Lsc     = 1.					# m
	psc		= Δρg * lc				# Pa
	tsc		= η_ϕbg_min / psc	# s
	println("lc : ",lc,";	tsc : ",tsc)
	# show_caracteristic
	η_r		= η0_r .* exp.(- θ_r .* (T.-T0))	# matrix viscosity (Pa.s)
	ηf_r	= η0_f .* exp.(- θ_f .* (T.-T0))	# fluid viscosity (Pa.s)
	η_ϕbg	= η_r / ϕ_bg					# Pa.s	
	Lc		= sqrt.(k_r ./ ηf_r .* η_r ./ ϕ_bg)	# m
	Tsc		= η_ϕbg ./ (Δρg .* Lc)	# s
	show_caracteristic(Lc, Tsc, xc, dir_plot,save_plots)
	# numerics time
	nt 		= 5e6
	dt_meca = 1e-6 * tsc 
	nvis	= nt/50  
	real_time = 0.0
	dt_diff	= 0.0
	dt_adv	= 0.0
	dt		= 0.0
	if do_p_vel dt_vel	= 0.0 end
	# init porosity
	ϕ		= @. ϕ_bg + ϕA * exp(-((xc + box_top + lx * loc_ϕ) / w) ^ 2)
	#ϕ[1] 	= 1e-7
	ϕ_init	= copy(ϕ)
	dϕdt    = zeros(nx - 2)
	ϕ_vel 	= Float64[] # ElasticArray{Float64}(undef, 1, 0)
	qD_vel	= Float64[]
	lim_bubb = Int(round(0.75 * nx))
	ϕ_heigh = bubble_heigh(ϕ, lim_bubb)
	# init pressure
	Pe		= zeros(nx)
	dPedt   = zeros(nx - 2)
	qD		= zeros(nx - 1)
	save_dt = zeros(Int(nt))
	# action
	for it = 1:Int(nt)
		# prorosity dependance
		ρCpT = ρf .* Cpf .* ϕ .* (1 .- ϕ) .+ ρs .* Cps
		ρCpf = (ρf .* Cpf)
		λT = λf .* avx(ϕ) .+ λs .* (1 .- avx(ϕ))
		χT = λT ./ avx(ρCpT)
		#test eta Frank
		η_r		= η0_r .* exp.(- θ_r .* (T.-T0))	# matrix viscosity (Pa.s)
	    ηf_r	= η0_f .* exp.(- θ_f .* (T.-T0))	# fluid viscosity (Pa.s)
		η_ϕbg	= η_r / ϕ_bg					# Pa.s	
		k_ηf	= k_r ./ ηf_r					# m3.s2.kg-1	
		Lc		= sqrt.(k_r ./ ηf_r .* η_r ./ ϕ_bg)	# m
		Tsc		= η_ϕbg ./ (Δρg .* Lc)	# s
		# Fluid flux
		qD   .= avx(k_ηf .* ϕ.^npow) .* (.- diff(Pe) ./ dx .- Δρg)
    	# time step 
		dt_adv		= 0.1*dx/maximum(abs.(qD))
		dt_diff		= 0.49*dx^2 / maximum(χT)
		dt = min(dt_adv,dt_meca,dt_diff)/10
		real_time   += dt 
		save_dt[it] = dt
		if do_p_vel dt_vel += dt end
		# Compute meca
		Compute_meca(ϕ, Pe, dPedt, qD, dϕdt, Δρg,k_ηf, η_r, β_f, npow, dx, dt, melt_input)
		# Compute themperature
		Compute_T(T, dU_dt, qD, ρCpT, λT, ρCpf, Grad_T, loc_T, T0, nx, dx, dt, dt_diff) 
		#Break ?
		do_break = break_conditions(ϕ,ϕ_init,xc,qD,T,T_init,real_time/year2sec,it,nx,lim_bubb,dx,box_top,year2sec)
		if do_break
			nt = it
			println("viscosity of the matrix : ",η_r[ϕ_heigh],";	viscosity of the melt : ", ηf_r[ϕ_heigh])
			ploting(ϕ,ϕ_init,T,T_init,ηf_r,η_r,xc,real_time,ϕ_heigh,it,year2sec,save_plots,do_p_bb_h,dir_plot)
			break
		end
		# visualisation 
		if (it % nvis) == 0 && do_visu
			old_ϕ_heigh = ϕ_heigh
			ϕ_heigh = bubble_heigh(ϕ, lim_bubb)
			if do_p_vel 
				push!(ϕ_vel, (ϕ_heigh - old_ϕ_heigh) * dx / (dt_vel / year2sec)) 
				push!(qD_vel, qD[ϕ_heigh])
				dt_vel = 0.0
			end		
			ploting(ϕ,ϕ_init,T,T_init,ηf_r,η_r,xc,real_time,ϕ_heigh,it,year2sec,save_plots,do_p_bb_h,dir_plot)
		end
	end
	
	if do_p_vel 
		display(plot(ϕ_vel, nvis:nvis:nt, xlims=(min(-0.005,minimum(ϕ_vel)*-0.9),maximum(ϕ_vel)*1.1)))
	end
	
	#Save data
	if do_save
		ϕ_heigh = bubble_heigh(ϕ, lim_bubb)
		T_diff = bubble_heat(T, T_init, ϕ_heigh,ϕ)
		to_save = [ϕ_heigh, T_diff, k_r, η0_r, θ_r, η0_f, θ_f, dt, Int(nt), dx, Int(nx), lc, Δρg, minimum(η_ϕbg), psc ,tsc]
		saving_data(to_save, f_result)
	end
end

# some boolean
do_visu		= true		# to make plot
do_p_bb_h	= true		# to visualize the bubble on plot
save_plots	= false		# to save plots
do_save  	= false
do_p_itr	= false 	# to print iterations

do_p_vel	= false		# to plot the velocity over time
melt_input	= false

# set rock viscosity parameters
η0_r	= 1e15			# [1e15 1e16] # viscosity at 750°C
θ_r 	= 0.038			# coefficient with temperature in exponential 
# set fluid viscosity parameters
η0_f	= 1e2			# [1e-3 1e-2 1e-1 1e0 1e1 1e2]			# viscosity at 750°C
θ_f 	= 0.03			# coefficient with temperature in exponential 

perm 	= 1e-12				# matrix permeabiblity

f_result = []
dir_plot = []

t_beg = Base.time()

dU_1D(η0_r, θ_r, η0_f, θ_f, perm, do_visu , do_p_bb_h, save_plots, do_save, do_p_itr, melt_input, dir_plot, f_result)

t_run = Base.time()-t_beg
if t_run>3600 Hours_run=round(t_run/3600,RoundDown) else Hours_run = 0 end
if t_run>60 Min_run=round(t_run/60-60*Hours_run,RoundDown) else Min_run = 0 end
secs_run=t_run-(60*Min_run+3600*Hours_run)
@printf("\n temps de calcule: %.0f h, %.0f min, %.0f sec \n", Hours_run, Min_run, secs_run)
