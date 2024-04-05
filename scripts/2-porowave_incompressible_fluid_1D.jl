using Plots, Plots.Measures, Printf
default(xmirror = true,size=(900, 600), framestyle=:box, label=false, grid=false, margin=10mm, lw=4, labelfontsize=20, tickfontsize=16, titlefontsize=20)

@views av(A) = 0.5 .* (A[1:end-1] .+ A[2:end])

@views function porowave_1D_fluid(do_visu)
lx   = 50        # longeur du modèle
kμ0  = 1         # permeabilite initial m2, viscosité fluide Pa s
ϕ0   = 0.01      # porosité du fond
npow = 3
Δϕ   = 0.1       # différence de porosité entre le fond et une bulle de porosité
ρfg  = 1         # masse volumique eau kg m-3, gravite terrestre m s-2
ρsg  = 2         # masse volumique rocks kg m-3, gravite terrestre m s-2
Δρg  = ρsg - ρfg # différence des masses volumiqueskg m-3, gravite terrestre m s-2
βϕ   = 0.1       # compressibilité du fluide 
η    = 1         # matrix bulk viscosity
# numerics
nx   = 200       # nombre de noeuds de modèles en x
dx   = lx / nx   # pas d'espace
xc   = LinRange((1:nx) .* dx)   # coordonnées en x des noeuds du modèles
nt   = 1e6       # nombre de pas de temps
dt   = 1e-3        # pas de temps
nvis =  nt / 100
# initialisation
ϕ      = @. Δϕ  * exp(-(xc - (lx*3/4 + dx / 2))^2) + ϕ0
Pe     = zeros(nx)
qDx    = zeros(nx - 1)
dPedt  = zeros(nx - 2)
dϕdt   = zeros(nx - 2)
#loop
for it = 1:nt 
    qDx   .= kμ0 .* av(ϕ).^npow .* (.- diff(Pe) ./ dx .+ Δρg)
    dϕdt  .= diff(qDx) ./ dx
    dPedt .= .- dϕdt .- ϕ[2:end-1] .* Pe[2:end-1] ./ η
    Pe[2:end-1] .+= dt ./ ϕ[2:end-1] ./ βϕ .* dPedt     
    ϕ[2:end-1]  .+= dt .* dϕdt
    # visualisation
    if (it % nvis) == 0 && do_visu
        p1 = plot(ϕ, xc, title="ϕ",yaxis= :flip, xlims=(0,0.4));
        p2 = plot(Pe, xc, title="Pe",yaxis= :flip,xlims=(-0.15,0.15));
        display(plot(p1,p2,plot_title=@sprintf(" time : %1.1e",dt*it/(365*24*3600))))
    end
end

end

porowave_1D_fluid(true)