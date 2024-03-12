function J = etc3(x_C3,h_m,qh2_m,q_m,pi_m,dG_H,dPsi_m,cox_i,cred_i,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electron Transport Chain Complex 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dG_C3op = - 35 + 2.*p.RT*log(h_m./1e-7)-p.RT.*log(qh2_m./q_m);

Pi3 =  ( 1 + pi_m ./ p.k_Pi3)./ (1 + pi_m ./ p.k_Pi4);

E_C3 = exp(-( dG_C3op + 4 .* dG_H - 2 .* p.F .* dPsi_m )./(2 .* p.RT ) );

J = x_C3 .* Pi3 .* (E_C3 .* cox_i - cred_i);
end