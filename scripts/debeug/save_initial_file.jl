
# some boolean
do_visu		= true		# to make plot
do_p_bb_h	= true		# to visualize the bubble on plot
save_plots	= false		# to save plots
do_p_vel	= false		# to plot the velocity over time
do_p_itr	= true 	# to print iterations
melt_input	= false
#set constant viscosity ?
cte_η_f = false			# to set a contant fluid viscosity
cte_η_r = false			# to set a contant rock viscosity
# set rock viscosity parameters
if cte_η_r
	η0_r 	= 1e22 		# constant rock viscosity 
	θ_r 	= 0.0		# do not change
else 
	η0_r 	= 1e15	# viscosity at 0°C
	θ_r 	= 0.038		# coefficient with temperature in exponential 
end 
# set fluid viscosity parameters
if cte_η_f
	η0_f	= 1e2		# constant fluid viscosity
	θ_f 	= 0.0		# do not change
else 
	η0_f	= 1e2		# viscosity at 0°C
	θ_f 	= 0.03		# coefficient with temperature in exponential 
end 

perm 	= 1e-11 		# matrix permeabiblity


t_beg = Base.time()

Temperature_E(do_visu, do_p_itr, do_p_bb_h, do_p_vel,save_plots,melt_input,cte_η_f,cte_η_r, perm, η0_r, θ_r, η0_f, θ_f)


t_run = Base.time()-t_beg
if t_run>3600 Hours_run=round(t_run/3600,RoundDown) else Hours_run = 0 end
if t_run>60 Min_run=round(t_run/60-60*Hours_run,RoundDown) else Min_run = 0 end
secs_run=t_run-(60*Min_run+3600*Hours_run)
@printf("\n temps de calcule: %.0f h, %.0f min, %.0f sec \n", Hours_run, Min_run, secs_run)
