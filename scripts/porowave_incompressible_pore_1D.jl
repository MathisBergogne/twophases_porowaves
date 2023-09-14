using Plots, Plots.Measures, Printf
default(xmirror = true,size=(1200, 800), framestyle=:box, label=false, grid=false, margin=10mm, lw=6, labelfontsize=20, tickfontsize=20, titlefontsize=24)

@views av(A) = 0.5 .* (A[1:end-1] .+ A[2:end])

@views function porowave_1D_pore(do_visu)
lx   = 50        # longeur du modèle
kμ0  = 1         # permeabilite initial m2, viscosité fluide Pa s
ϕ0   = 0.01      # porosité du fond
npow = 3
Δϕ   = 0.1       # différence de porosité entre le fond et une bulle de porosité
ρfg  = 1         # masse volumique eau kg m-3, gravite terrestre m s-2
ρsg  = 2         # masse volumique rocks kg m-3, gravite terrestre m s-2
Δρg  = ρsg - ρfg # différence des masses volumiqueskg m-3, gravite terrestre m s-2
Bf   = 0.1       # compressibilité du fluide 
η    = 1         # matrix bulk viscosity
# numerics
nx   = 200       # nombre de noeuds de modèles en x
dx   = lx / nx   # pas d'espace
xc   = LinRange((1:nx) .* dx)   # coordonnées en x des noeuds du modèles
nt   = 5e6       # nombre de pas de temps
dt   = 1e-3        # pas de temps
nvis = nt / 100
maxiter = 20nx
ncheck  = ceil(Int,0.05nx)
ϵtol    = 1e-8
# initialisation
ϕ      = @. Δϕ  * exp(-(xc - (lx*3/4 + dx / 2))^2) + ϕ0
Pe     = zeros(nx)
qDx    = zeros(nx - 1)
dPedt  = zeros(nx - 2)
dϕdt   = zeros(nx - 2)
# visualisation
p1 = plot(ϕ, xc, title="ϕ",yaxis= :flip, xlims=(0,0.4))
p2 = plot(Pe, xc, title="Pe",yaxis= :flip,xlims=(-0.15,0.15)) 
display(plot(p1,p2))

#loop
for it = 1:nt 
    qDx   .= kμ0 .* av(ϕ).^npow .* (.- diff(Pe) ./ dx .+ Δρg)
    dϕdt  .= .- ϕ[2:end-1] .* Pe[2:end-1] ./ η
    dPedt .= .- diff(qDx) ./ dx .- ϕ[2:end-1] .* Pe[2:end-1] ./ η
    Pe[2:end-1] .+= dt ./ ϕ[2:end-1] ./ Bf .* dPedt
    ϕ[2:end-1]  .+= dt .* dϕdt
    if isnan(ϕ[2]) break end
    # visualisation
    if (it % nvis) == 0 && do_visu
        p1 = plot(ϕ, xc, title="ϕ",yaxis= :flip, xlims=(0,0.4))
        p2 = plot(Pe, xc, title="Pe",yaxis= :flip,xlims=(-0.15,0.15)) 
        display(plot(p1,p2))
    end
end

end

porowave_1D_pore(true)
