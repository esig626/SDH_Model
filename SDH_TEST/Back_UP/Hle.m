function J = Hle(gHLe,dPsi_m,h_m,p)

J = gHLe * dPsi_m * (p.H_i*exp(p.F*dPsi_m/p.RT)-h_m)/(exp(p.F*dPsi_m/p.RT)-1);

end