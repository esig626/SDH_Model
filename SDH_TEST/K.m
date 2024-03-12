function J = K(x_K,p,dPsi_m,k_m)

J = x_K * dPsi_m * (p.K_i*exp(p.F*dPsi/p.RT)-k_m)/(exp(p.F*dPsi_m/p.RT)-1);


end