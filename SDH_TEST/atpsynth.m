function J= atpsynth(x_F1,pi_m, adp_m, atp_m, dG_H,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ATP-Synthase (Complex V)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_ATP = exp(-( p.dG_F1o - p.n_A .* dG_H) ./ p.RT );

E_5 = pi_m .* (p.K_DD ./ p.K_DT);

J = x_F1 .* (E_ATP .* E_5 .* adp_m - atp_m); 

end