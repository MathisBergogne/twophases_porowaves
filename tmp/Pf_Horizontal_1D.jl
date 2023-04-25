using Plots, Plots.Measures
default(size=(1200, 800), framestyle=:box, label=false, grid=false, margin=10mm, lw=6, labelfontsize=20, tickfontsize=20, titlefontsize=24)

# physics
lx   = 10        # longeur du modèle
kμf0 = 1         # permeabilite initial m2
ϕ    = 0.05      # porosité
# numerics
nx   = 100       # nombre de noeuds de modèles en x
dx   = lx / nx   # pas d'espace
xc   = LinRange((1:nx) .* dx)   # coordonnées en x des noeuds du modèles
nt   = 1e6       # nombre de pas de temps
dt   = dx^2 / kμf0 / 2.1        # pas de temps
# initialisation
qx   = zeros(nx - 1)
Pf   = @. 1.0 .* exp(-(xc - (lx / 2 + dx / 2))^2)
Pf_ini = copy(Pf)
#loop
for it = 1:nt
    qx .= .-kμf0 .* (diff(Pf) ./ dx)
    Pf[2:end-1] .-= dt .* diff(qx) ./ dx
    Pf[1] = Pf[2]
    Pf[end] = Pf[end-1]
end
# visualisation
plot(xc,[Pf,Pf_ini])

