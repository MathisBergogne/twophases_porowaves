using Plots, Plots.Measures, Printf
default(size=(1200, 800))#, xmirror = true, framestyle=:box, label=false, grid=false, margin=10mm, lw=6, labelfontsize=20, tickfontsize=20, titlefontsize=24)

@views avx(A) = 0.5 .* (A[1:end-1] .+ A[2:end])

@views function porowave_1D_speudotime_auto(do_visu)
lx   = 100        # longeur du modèle
ϕ0   = 0.01      # porosité du fond
npow = 3         # puissnace m pour ϕ
Δϕ   = 0.05       # différence de porosité entre le fond et une bulle de porosité
ρfg  = 1         # masse volumique eau kg m-3, gravite terrestre m s-2
ρsg  = 2         # masse volumique rocks kg m-3, gravite terrestre m s-2
Δρg  = ρsg - ρfg # différence des masses volumiqueskg m-3, gravite terrestre m s-2
β    = 1e-10     # compressibilité du fluide 
η    = 1         # matrix bulk viscosity
kμ0  = 1         # permeabilite initial m2, viscosité fluide Pa s
# numerics
nx   = 1000      # nombre de noeuds de modèles en x
dx   = lx / nx    # pas d'espace
xc   = LinRange(-dx/2, lx+dx/2,nx)   # coordonnées en x des noeuds du modèles
nt   = 1e2#6       # nombre de pas de temps
dt   = 1e-3        # pas de temps
nvis = nt / 100
maxiter = 20nx
ncheck  = ceil(Int,0.05nx)
ϵtol    = 1e-8 
cfl     = 1 / 1.1
# initialisation
ϕ      = @. Δϕ  * exp(-(xc - (lx*9/10 + dx / 2))^2 / 10.0) + ϕ0
Pe     = zeros(nx)
Pe_old = zeros(nx)
qD     = zeros(nx - 1)
RPe    = zeros(nx - 2)
RqD    = zeros(nx - 1)
η_ϕ    = ones(nx)
k_ηf   = zeros(nx)
lc_loc = zeros(nx)
re     = zeros(nx)
ρ_dτ   = zeros(nx)
vsdτ   = zeros(nx)
βf_dτ  = zeros(nx)
η_ϕτ   = zeros(nx)
k_ηfτ  = zeros(nx - 1)

# visualisation
p1 = plot(ϕ, xc, title="ϕ",yaxis= :flip,xlims=(0,0.11))
p2 = plot(Pe, xc, title="Pe",yaxis= :flip,xlims=(-4.5,4.5))
display(plot(p1,p2))

#loop
for it = 1:nt 
    Pe_old .= Pe
    iter = 1; err = 2ϵtol;
    while err >= ϵtol && iter <= maxiter
        η_ϕ    .= (ϕ0 ./ ϕ)
        k_ηf   .= (ϕ ./ ϕ0) .^ npow

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
        if iter % ncheck == 0
            err_qD = maximum(abs.(RqD))
            err_Pe = maximum(abs.(RPe))
            err    = max(err_Pe,err_qD)
        end
        iter += 1
    end
    #porosity update
    ϕ[2:end-1]  .+= dt .* Pe[2:end-1] ./ η_ϕ[2:end-1] 
    if isnan(ϕ[2]) break end
    # visualisation 
    if (it % nvis) == 0 && do_visu
        p1 = plot(ϕ, xc, title="ϕ",yaxis= :flip,xlims=(0,0.11))
        p2 = plot(Pe, xc, title="Pe",yaxis= :flip,xlims=(-4.5,4.5))
        display(plot(p1,p2))
    end
end

end

porowave_1D_speudotime_auto(true)
