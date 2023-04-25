using Plots,Plots.Measures
default(yflip = true,size=(1200,800),framestyle=:box,label=false,grid=false,margin=10mm,lw=6,labelfontsize=20,tickfontsize=20,titlefontsize=24)


#props
k0  = 1;     #permeabilite initial m2
phi = 0.05;      # porosité
mu  = 1;       #Pa s
rho = 1#2800      #masse volumique kg m-3
g   = 1#9.81      #gravite terrestre m s-2
#grid init 
lx  = 10;         #longeur du modèle
nx  = 100;        #nombre de noeuds de modèles en x
dx  = lx/nx;      #pas d'espace
xc  = LinRange(((1:nx).*dx));    #coordonnées en x des noeuds du modèles
nt  = nx^2/5;         #nombre de pas de temps         
dt  = dx^2/(k0*mu)/2.1;           #pas de temps
#array init
Qd = zeros(nx-1);
Pf   = @. 1.0 .* exp(-(xc - (lx / 2 + dx / 2))^2)
Pfcopy=copy(Pf);
#loop
for it = 1:nt
    Qd[:]  .= .-(k0/mu).*((diff(Pf))./dx .-rho*g)
    Pf[2:end-1] .-= dt.*diff(Qd)./dx
    Pf[1] = Pf[2] - rho*g *dx; 
    Pf[end] = rho*g *xc[end];
end

plot([Pf,Pfcopy],xc)

