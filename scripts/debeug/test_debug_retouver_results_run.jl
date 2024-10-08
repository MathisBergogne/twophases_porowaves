using Plots, ElasticArrays
using CSV, DataFrames, Printf

include("test_debug_retouver_results.jl")

# some boolean
do_visu		= true		# to make plot
do_p_bb_h	= true		# to visualize the bubble on plot
save_plots	= false		# to save plots
do_save  	= false
do_p_itr	= false 	    # to print iterations

do_p_vel	= false		# to plot the velocity over time
melt_input	= false

# set rock viscosity parameter 
η0_r 	= 1e15		# viscosity at 0°C
θ_r 	= 0.038		# coefficient with temperature in exponential 
# set fluid viscosity parameters
η0_f	= 1e-1		# viscosity at 0°C
θ_f 	= 0.03		# coefficient with temperature in exponential 

perm 	= 1e-12 		# matrix permeabiblity


@views function main(do_visu, do_p_bb_h, save_plots, do_save, do_p_itr, melt_input, η0_r, θ_r, η0_f, θ_f, perm)
	f_result = []
	dir_plot = []
	test_debug_refound(do_visu,do_p_itr,do_p_bb_h,do_p_vel,save_plots,do_save,melt_input,perm,η0_r,θ_r,η0_f,θ_f,dir_plot,f_result)
    
end

t_beg = Base.time()

main(do_visu, do_p_bb_h, save_plots, do_save, do_p_itr, melt_input, η0_r, θ_r, η0_f, θ_f, perm)

t_run = Base.time()-t_beg
if t_run>3600 Hours_run=round(t_run/3600,RoundDown) else Hours_run = 0 end
if t_run>60 Min_run=round(t_run/60-60*Hours_run,RoundDown) else Min_run = 0 end
secs_run=t_run-(60*Min_run+3600*Hours_run)
@printf("\n temps de calcule: %.0f h, %.0f min, %.0f sec \n", Hours_run, Min_run, secs_run)
