function value=qVT(rho,p,Ev)
%% ignore VTO for now (11/10/09)
global d kB x m
n     = rho/m.g;
T     = p./(n*kB);
nN2   = x.N2*n;
EvTg  = d.Ev./(exp(d.Ev./(kB*T))-1); % _J, equilibrium vibrational energy at T.g
value = nN2.*(Ev-EvTg)./tauVT(rho,p); % _J_m^-3/_s
end