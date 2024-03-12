function J = etc1(x_C1,nadh_m,nad_m,dG_H,q_m,qh2_m,h_m,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electron Transport Chain Complex 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dG_C1op = - 74 - p.RT*log(h_m/1e-7) - p.RT*log(q_m/qh2_m);

E_C1 = exp( - ( dG_C1op + 4 * dG_H)/p.RT);

J = x_C1 * ( E_C1 * nadh_m - nad_m);
end