using Plots, Plots.Measures
default(xmirror = true,size=(1200, 800), framestyle=:box, label=false, grid=false, margin=10mm, lw=3, labelfontsize=20, tickfontsize=15, titlefontsize=20)

@views function hydrostatic_1D()
# physics
lx   = 10        # longeur du modèle
kμf0 = 1         # permeabilite initial m2
ϕ    = 0.05      # porosité
ρfg  = 1         # masse volumique kg m-3, gravite terrestre m s-2
# numerics
nx   = 100       # nombre de noeuds de modèles en x
dx   = lx / nx   # pas d'espace
xc   = LinRange((1:nx) .* dx)   # coordonnées en x des noeuds du modèles
nt   = 1e5       # nombre de pas de temps
dt   = dx^2 / kμf0 / 2.1        # pas de temps
# initialisation
qx   = zeros(nx - 1)
Pf   = @. 1.0 .* exp(-(xc - (lx / 2 + dx / 2))^2)
# time loop
for it = 1:nt
    qx .= .-kμf0 .* (diff(Pf) ./ dx .- ρfg)
    Pf[2:end-1] .-= dt .* diff(qx) ./ dx
    Pf[1] = Pf[2] - ρfg * dx # conditions au bord no flux
    Pf[end] = ρfg*xc[end] # fix pressure
end
# visualisation
p1 = plot(Pf, xc, title="Pf",yaxis= :flip)
p2 = plot(qx, xc[2:end], title="qx")
plot(p1, p2)

end

hydrostatic_1D()