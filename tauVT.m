function value=tauVT(rho,p)
global kB m
n     = rho/m.g;
T     = p./(n*kB);
value = 6.5e-4*exp(137./T.^(1/3))./p; % _s, [Mnatskanyan:1985]
end