using Plots, Plots.Measures, Printf
default(size=(1200, 800), framestyle=:box , label=false,  margin=5mm , lw=2, labelfontsize=20, tickfontsize=15, titlefontsize=18)

@views avx(A) = 0.5 .* (A[1:end-1] .+ A[2:end])

@views function porowave_1D_speudotime_auto(do_visu)
lx   = 100        # longeur du modèle
ϕ0   = 0.01      # porosité du fond
npow = 3         # puissnace m pour ϕ
Δϕ   = 0.1       # différence de porosité entre le fond et une bulle de porosité
ρfg  = 1         # masse volumique eau kg m-3, gravite terrestre m s-2
ρsg  = 2         # masse volumique rocks kg m-3, gravite terrestre m s-2
Δρg  = ρsg - ρfg # différence des masses volumiqueskg m-3, gravite terrestre m s-2
β    = 1e-10     # compressibilité du fluide 
η    = 1         # matrix bulk viscosity
kμ0  = 1         # permeabilite initial m2, viscosité fluide Pa s
ΔT   = 200.0
T0   = 500.0
Ra   = 1000.0
αρgx,αρgy = 0.0,1.0
αρg       = sqrt(αρgx^2+αρgy^2)
# numerics
nx   = 100      # nombre de noeuds de modèles en x
dx   = lx / nx    # pas d'espace
xc   = LinRange(-dx/2, lx+dx/2,nx)   # coordonnées en x des noeuds du modèles
nt   = 1e2       # nombre de pas de temps
dt   = 1e-3        # pas de temps
dτ   = 0.01
nvis = nt / 100
maxiter = 20nx
ncheck  = ceil(Int,0.05nx)
ϵtol    = 1e-8
re_D    = 4π
cfl     = 1 / 1.1
# initialisation
ϕ      = @. Δϕ  * exp(-(xc - (lx*8/10 + dx / 2))^2 / 10.0) + ϕ0
φ_init = copy(ϕ)
T      = @. T0 + ΔT * exp(-(xc - (lx*8/10 + dx / 2))^2 / 5.0)
T_init = copy(T)
T_old  = zeros(nx)
r_T    = zeros(nx - 2)
Pe     = zeros(nx)
Pe_old = zeros(nx)
qD     = zeros(nx - 1)
RPe    = zeros(nx - 2)
RqD    = zeros(nx - 1)
dTdt   = zeros(nx - 2)
qT     = zeros(nx - 1)
η_ϕ    = zeros(nx)
k_ηf   = zeros(nx)
lc_loc = zeros(nx)
re     = zeros(nx)
ρ_dτ   = zeros(nx)
vsdτ   = zeros(nx)
βf_dτ  = zeros(nx)
η_ϕτ   = zeros(nx)
k_ηfτ  = zeros(nx - 1)
λ_ρCp  = zeros(nx)
re_T   = zeros(nx)
θ_dτ_T = zeros(nx)
β_dτ_T = zeros(nx)
Count_iter = Float64[]
#loop
for it = 1:nt 
    T_old  .= T
    Pe_old .= Pe
    iter = 1; err = 2ϵtol;iter_evo = Float64[]; err_evo = []
    while err >= ϵtol && iter <= maxiter
        η_ϕ    .= (ϕ0 ./ ϕ)
        k_ηf   .= (ϕ ./ ϕ0) .^ npow

        λ_ρCp   = 1 ./ Ra .* (αρg .* k_ηf .* ΔT .* lx ./ ϕ)
        re_T    = π .+ sqrt.(π^2 .+ lx^2 ./ λ_ρCp ./ dt)
        θ_dτ_T  = lx ./ re_T ./ cfl ./ dx
        β_dτ_T  = (re_T .* λ_ρCp) ./ (cfl .*dx .* lx)

        lc_loc .= sqrt.(k_ηf .* η_ϕ)
        re     .= π .+ sqrt.(π .^ 2 .+ (lx ./ lc_loc) .^ 2)
        ρ_dτ   .= re .* 1.0 ./ (1.0 ./ η_ϕ .+ β ./ dt) ./ lx ./ dx ./ cfl
        vsdτ    = cfl * dx
        βf_dτ  .= 1 / vsdτ^2 ./ ρ_dτ
        η_ϕτ   .= 1.0 ./ (βf_dτ .+ 1.0 ./ η_ϕ .+ β ./ dt)
        k_ηfτ  .= 1.0 ./ (avx(ρ_dτ) .+ 1.0 ./ avx(k_ηf))

        # fluid pressure update
        RPe   .= .+ diff(qD) ./ dx .- Pe[2:end-1] ./ η_ϕ[2:end-1] .- β .* (Pe[2:end-1] .- Pe_old[2:end-1]) ./ dt 
        Pe[2:end-1] .+= η_ϕτ[2:end-1] .* RPe
        RqD     .= .-qD ./ avx(k_ηf) .+ diff(Pe) ./ dx .+ Δρg
        qD     .+= k_ηfτ .* RqD
        # temperature update
        qT  .-= (qT .+ avx(λ_ρCp) .*(diff(T)./dx))./(1.0 .+ avx(θ_dτ_T))
        dTdt[2:end-1] .= (T[3:end-2] .- T_old[3:end-2])./dt #.+ 
            (max.(qD[2:end-2],0.0).*diff(T[2:end-2])./dx .+
             min.(qD[3:end-1],0.0).*diff(T[3:end-1])./dx )
        T[2:end-1] .-= (dTdt .+ diff(qT)./dx)./(1.0/dt .+ β_dτ_T[2:end-1])
        #T[[1,end]]  .= T[[2,end-1]]

        if iter % ncheck == 0
            err_qD = maximum(abs.(RqD))
            err_Pe = maximum(abs.(RPe))
            r_T   .= dTdt .+ diff(qT)./dx
            err_T  = maximum(abs.(r_T))
            err    = max(err_Pe,err_T,err_qD)
            push!(err_evo,err)
            push!(iter_evo, iter / nx)
        end
        iter += 1
    end
    push!(Count_iter, iter / nx)
    #porosity update
    ϕ[2:end-1]  .+= dt .* Pe[2:end-1] ./ η_ϕ[2:end-1] 
    if isnan(ϕ[2]) break end
    # visualisation 
    if (it % nvis) == 0 && do_visu
        p1 = plot([ϕ,φ_init], xc, title="ϕ",yaxis= :flip, xlims=(0,0.11), label=["ϕ", "φ_init"])
        p2 = plot(Pe, xc, title="Pe",yaxis= :flip,xlims=(-4.5,4.5)) 
        p3 = plot([T,T_init], xc, title="T",yaxis= :flip, xlims=(490,710))
        p4 = plot(iter_evo,err_evo,title="iter / nx",yaxis=:log10, minorgrid=true, marker=:circle)
        p5 = plot(1:it,Count_iter,title="iter / nx", marker=:circle)
        display(plot(p1,p2,p3,p4,p5))
    end
end

end

porowave_1D_speudotime_auto(true)
