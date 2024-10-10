using Plots, ElasticArrays
using CSV, DataFrames, Printf

include("4-Porowave_1D_dU.jl")

# some boolean
do_visu		= true		# to make plot
do_p_bb_h	= true		# to visualize the bubble on plot
save_plots	= false		# to save plots
do_save  	= false
do_p_itr	= false 	# to print iterations
do_p_var	= true

do_p_vel	= false		# to plot the velocity over time
melt_input	= false

# set rock viscosity parameters
	η0_R	= 1e15 # [1e15 1e16 1e17 1e18 1e19]	# [1e15 1e16] # viscosity at 750°C
	θ_r 	= 0.038			# coefficient with temperature in exponential 
# set fluid viscosity parameters
	η0_F	= 1e2 # [1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3] # [1e-3 1e-2 1e-1 1e0 1e1 1e2]			# viscosity at 750°C
	θ_f 	= 0.03			# coefficient with temperature in exponential 

Perm 	= [1e-11 1e-12 1e-13]  			# matrix permeabiblity


@views function main(do_visu, do_p_bb_h, save_plots, do_save, do_p_itr, do_p_var, do_p_vel, melt_input, η0_R, θ_r, η0_F, θ_f, Perm)
	cd("/home/bergogne/Documents/6-Julia/twophases_porowaves/etud_param")
	name = @sprintf("correction_dt_no_melt-input/%1.0e-%1.0e-Vf_%1.0e-%1.0e-Vr",minimum(η0_F),maximum(η0_F),minimum(η0_R),maximum(η0_R))
	if do_save
		f_result = name * ".csv" # "/etud_param_Porowaves_1D_dU.csv"
		results_file(f_result, do_save)
	end 
	
	for i in eachindex(η0_R)
		η0_r = η0_R[i]
		if do_p_var @printf("\nη0_r = %1.2e; \n", η0_r) end
		for j in eachindex(η0_F)
			η0_f = η0_F[j]
			if do_p_var @printf("  η0_f = %1.2e; \n", η0_f) end
			for k in eachindex(Perm)
				perm = Perm[k]
				if do_p_var @printf("    perm = %1.2e; \n", perm) end
				@printf("		η0_r = %1.2e;	η0_f = %1.2e; \n", η0_r, η0_f)
				if save_plots
					dir_plot = @sprintf("correction_dt_no_melt-input/Runs_%1.0ek_%1.0efv_%1.0esv",perm,η0_f,η0_r)
					if isdir(dir_plot)
						save_plots = false
						dir_plot = []
					else
						mkdir(dir_plot)
					end
				else
					dir_plot = []
				end 
				dU_1D(η0_r, θ_r, η0_f, θ_f, perm, do_visu , do_p_bb_h, save_plots, do_save, do_p_itr, melt_input, dir_plot, f_result)
			end
		end
	end
end

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

t_beg = Base.time()

main(do_visu, do_p_bb_h, save_plots, do_save, do_p_itr, do_p_var, do_p_vel, melt_input, η0_R, θ_r, η0_F, θ_f, Perm)

t_run = Base.time()-t_beg
if t_run>3600 Hours_run=round(t_run/3600,RoundDown) else Hours_run = 0 end
if t_run>60 Min_run=round(t_run/60-60*Hours_run,RoundDown) else Min_run = 0 end
secs_run=t_run-(60*Min_run+3600*Hours_run)
@printf("\n temps de calcule: %.0f h, %.0f min, %.0f sec \n", Hours_run, Min_run, secs_run)
