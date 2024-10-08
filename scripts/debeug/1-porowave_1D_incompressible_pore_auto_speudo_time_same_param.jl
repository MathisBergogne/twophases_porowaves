using Plots, Plots.Measures, Printf
default(size=(900, 600), framestyle=:box , label=false,  margin=5mm , lw=2, labelfontsize=20, tickfontsize=10, titlefontsize=16)

@views avx(A) = 0.5 .* (A[1:end-1] .+ A[2:end])

@views function porowave_1D_speudotime_auto(do_visu)
auto_speudo_time = true
lx   = 10000		# longeur du modèle
ϕ0   = 0.01			# porosité du fond
npow = 3			# puissnace m pour ϕ
Δϕ   = 0.1			# différence de porosité entre le fond et une bulle de porosité
ρfg  = 2000*10		# masse volumique eau kg m-3, gravite terrestre m s-2
ρsg  = 3000*10 		# masse volumique rocks kg m-3, gravite terrestre m s-2
Δρg  = ρsg - ρfg	# différence des masses volumiqueskg m-3, gravite terrestre m s-2
η    = 1e15		    # matrix bulk viscosity
kμ0  = 1e-12/1e3			# permeabilite initial m2, viscosité fluide Pa s
if auto_speudo_time
    β    = 1e-10
    println("lc : ",sqrt(kμ0*η/ϕ0^2))
else 
    β    = 1e-10
    println("lc : ",sqrt(kμ0*η/ϕ0^2))
end		    # compressibilité du fluide 
# numerics
nx   = 1000			# nombre de noeuds de modèles en x
dx   = lx / nx		# pas d'espace
xc   = LinRange(-dx/2, lx+dx/2,nx)   # coordonnées en x des noeuds du modèles
nt   = 1e3			# nombre de pas de temps
dt   = 1e9			# pas de temps
dτ   = 1e5
nvis = nt / 100
maxiter = 20nx
ncheck  = ceil(Int,0.05nx)
etol    = 1e-8 
cfl     = 1 / 1.1
# initialisation
ϕ      = @. Δϕ  * exp(-(xc - (lx*9/10 + dx / 2))^2 / 10000.0) + ϕ0 # exp(-(xc - (lx*9/10 + dx / 2))^2 / 10.0) + ϕ0
φ_init = copy(ϕ) 
φ_old = copy(ϕ)   
Pe     = zeros(nx)
Pe_old = zeros(nx)
η_ϕ    = ones(nx)

###
qD     = zeros(nx - 1)
RPe    = zeros(nx - 2)
RqD    = zeros(nx - 1)
k_ηf   = zeros(nx)
lc_loc = zeros(nx)
re     = zeros(nx)
ρ_dτ   = zeros(nx)
vsdτ   = zeros(nx)
βf_dτ  = zeros(nx)
η_ϕτ   = zeros(nx)
k_ηfτ  = zeros(nx - 1)
###
qD     = zeros(nx - 1)
dPedτ  = zeros(nx - 2)
dϕdτ   = zeros(nx - 2)

# visualisation
p1 = plot(ϕ, xc, title="ϕ",yaxis= :flip,xlims=(0,0.11))
p2 = plot(Pe, xc, title="Pe",yaxis= :flip,xlims=(-4.5,4.5))
display(plot(p1,p2))

#loop
for it = 1:nt 
    Pe_old .= Pe
    if auto_speudo_time
        meca_speudotime_auto(ϕ, qD, Pe, Pe_old, RPe, RqD, kμ0, η, η_ϕ, k_ηf, ϕ0, Δρg, npow, lc_loc, re, ρ_dτ, vsdτ, βf_dτ, η_ϕτ, k_ηfτ, β, lx, dx, dt, cfl, etol, maxiter, ncheck)
    else
        meca_speudotime(ϕ, φ_old, Pe, dPedτ, qD, dϕdτ, Δρg, kμ0 , η, β, npow, dx, nx, dt, dτ, etol, maxiter, ncheck)
    end
    #porosity update
    #ϕ[end] 		= 1e-4
    if isnan(ϕ[2]) 
        println("break")
        break 
    end
    # visualisation 
    if (it % nvis) == 0 && do_visu
        p1 = plot([ϕ,φ_init], xc, title="ϕ",yaxis= :flip,xlims=(0,0.30))
        p2 = plot(Pe, xc, title="Pe",yaxis= :flip)
        p3 = plot(qD, avx(xc), title="qD",yaxis= :flip)
        display(plot(p1,p2,p3,plot_title=@sprintf(" time : %1.1e y",it*dt/(365*24*3600)); layout=(1, 3)))
        #println("max/min qD : ",maximum(qD)," / ",minimum(qD))
    end
end

end

@views function meca_speudotime_auto(ϕ, qD, Pe, Pe_old, RPe, RqD, kμ0, η, η_ϕ, k_ηf, ϕ0, Δρg, npow, lc_loc, re, ρ_dτ, vsdτ, βf_dτ, η_ϕτ, k_ηfτ, β, lx, dx, dt, cfl, etol, maxiter, ncheck)
    iter = 1; err = 2etol;
    while err >= etol && iter <= maxiter
        η_ϕ    .= η ./ϕ0 .* (ϕ0 ./ ϕ)
        k_ηf   .= kμ0 .* (ϕ ./ ϕ0) .^ npow

        lc_loc .= sqrt.(k_ηf .* η_ϕ)
        re     .= π .+ sqrt.(π .^ 2 .+ (lx ./ lc_loc) .^ 2)
        ρ_dτ   .= re .* 1.0 ./ (1.0 ./ η_ϕ .+ β ./ dt) ./ lx ./ dx ./ cfl
        vsdτ    = cfl * dx
        βf_dτ  .= 1 / vsdτ^2 ./ ρ_dτ
        η_ϕτ   .= 1.0 ./ (βf_dτ .+ 1.0 ./ η_ϕ .+ β ./ dt)
        k_ηfτ  .= 1.0 ./ (avx(ρ_dτ) .+ 1.0 ./ avx(k_ηf))

        # fluid pressure update
        RqD     .= .-qD ./ avx(k_ηf) .+ diff(Pe) ./ dx .+ Δρg
        qD     .+= k_ηfτ .* RqD
        RPe   .= .+ diff(qD) ./ dx .- Pe[2:end-1] ./ η_ϕ[2:end-1] .- β .* (Pe[2:end-1] .- Pe_old[2:end-1]) ./ dt 
        Pe[2:end-1] .+= η_ϕτ[2:end-1] .* RPe
        # temperature update
        if iter % ncheck == 0
            err_qD = maximum(abs.(RqD))
            err_Pe = maximum(abs.(RPe))
            err    = max(err_Pe,err_qD)
        end
        iter += 1
    end
    ϕ[2:end-1]  .+= dt .* Pe[2:end-1] ./ η_ϕ[2:end-1] 
end

@views function meca_speudotime(ϕ, φ_old, Pe, dPedτ, qD, dϕdτ, Δρg, kμ0 , η, Bf, npow, dx, nx, dt, dτ, etol, maxiter, ncheck)
    φ_old = ϕ
    iter = 1;err = 2etol;iter_evo = Float64[]; err_evo = []
    while err >= etol && iter <= maxiter
        qD    .= kμ0 .* avx(ϕ).^npow .* (.- diff(Pe) ./ dx .+ Δρg)
        dϕdτ  .= .- ϕ[2:end-1] .* Pe[2:end-1] ./ η
        dPedτ .= .- diff(qD) ./ dx .+ dϕdτ # .- ϕ[2:end-1] .* Pe[2:end-1] ./ η
        Pe[2:end-1] .+= dτ ./ ϕ[2:end-1] ./ Bf .* dPedτ
        ϕ[2:end-1]  .+= dτ .* dϕdτ
        if iter % ncheck == 0
            err = maximum(abs.(diff(Bf.*diff(ϕ)./dx)./dx .- (ϕ[2:end-1] .- φ_old[2:end-1])./dt))
            push!(iter_evo,iter/nx); push!(err_evo,err)
        end
        iter += 1
    end
end

porowave_1D_speudotime_auto(true)
