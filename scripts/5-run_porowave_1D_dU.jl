using Plots, ElasticArrays
using CSV, DataFrames, Printf

include("5-Porowave_1D_dU_auto_speudo_time.jl")

# some boolean
do_visu		= true		# to make plot
do_p_bb_h	= true		# to visualize the bubble on plot
save_plots	= false		# to save plots
do_save  	= false
do_p_itr	= false 	# to print iterations

do_p_vel	= false		# to plot the velocity over time
melt_input	= false

# set rock viscosity parameters
	η0_R	= 5e17			# [1e15 1e16] # viscosity at 750°C
	θ_r 	= 0.04			# coefficient with temperature in exponential 
# set fluid viscosity parameters
	η0_F	= 1e3			# [1e-3 1e-2 1e-1 1e0 1e1 1e2]			# viscosity at 750°C
	θ_f 	= 0.03			# coefficient with temperature in exponential 

Perm 	= 1e-13				# matrix permeabiblity


@views function main(do_visu, do_p_bb_h, save_plots, do_save, do_p_itr, melt_input, η0_R, θ_r, η0_F, θ_f, Perm)
	f_result = []
	if save_plots
		dir_plot = @sprintf("ref_visco-r_%1.1e_visco-f_%1.1e_perm_%1.1e",η0_R,η0_F,Perm) # []
	else
		dir_plot = []
	end
	if dir_plot != [] && save_plots
		mkdir(dir_plot)
	end	
	η0_r = η0_R
	η0_f = η0_F
	perm = Perm
	dU_1D(η0_r, θ_r, η0_f, θ_f, perm, do_visu , do_p_bb_h, save_plots, do_save, do_p_itr, melt_input, dir_plot, f_result)
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

main(do_visu, do_p_bb_h, save_plots, do_save, do_p_itr, melt_input, η0_R, θ_r, η0_F, θ_f, Perm)

t_run = Base.time()-t_beg
if t_run>3600 Hours_run=round(t_run/3600,RoundDown) else Hours_run = 0 end
if t_run>60 Min_run=round(t_run/60-60*Hours_run,RoundDown) else Min_run = 0 end
secs_run=t_run-(60*Min_run+3600*Hours_run)
@printf("\n temps de calcule: %.0f h, %.0f min, %.0f sec \n", Hours_run, Min_run, secs_run)
