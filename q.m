function value=q(r)
global q0 rs
value = q0*exp(-r.^2/rs^2);
end