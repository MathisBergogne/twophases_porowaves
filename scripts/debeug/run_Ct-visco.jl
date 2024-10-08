using Plots, Plots.Measures, Printf, ElasticArrays

include("../porowave_1D_with_T_from_U.jl")

t_tic = Base.time()

do_visu   = true
do_p_itr  = false 
do_save	  = false

perm 	= 1e-13 #[1e-12 2.5e-12 5e-12 7.5e-12 1e-11 1e-11 ];
visco_f = 1e5 #[1e4   1e4     1e4   1e4     1e4   4e3   ];
visco_r = 1e19 #[1e16  2.5e16  5e16  7.5e16  1e17  2.5e17];

lc		= sqrt(perm / visco_f * visco_r / 0.01)
println("lc : ", lc)

for i in eachindex(perm)
	if do_save
		dir_plots = @sprintf("../run_T/%1.0ek_%1.0efv_%1.0esv/",perm[i],visco_f[i],visco_r[i])
		if isdir(dir_plots)
			do_saving = false
		else
			mkdir(dir_plots)
			do_saving = do_save
		end
	else 
		dir_plots = []
		do_saving = do_save
	end
	Temperature_E(do_visu, do_p_itr, do_saving, perm[i], visco_f[i], visco_r[i],dir_plots)
end

t_rest = ((Base.time()-t_tic))
Hours=round(t_rest/3600,RoundDown)
Min=round(t_rest/60-60*Hours,RoundDown)
secs=t_rest-(60*Min+3600*Hours)
@printf("time : %.0f h, %.0f min, %.0f sec \n", Hours, Min, secs);
