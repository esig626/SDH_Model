function J = etc4(x_C4,cred_i,cox_i,dPsi_m,dG_H,o2,h_m,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electron Transfer Chain (Complex IV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dG_C4op = p.dG_C4o - 2 * p.RT * log(h_m / 1e-7)- p.RT/2*log(o2/1);

OX_4 =  1 / (1 + p.k_O2/o2 ) * cred_i / p.Ctot;

E_4 = exp(-( dG_C4op + 2 * dG_H)/(2 * p.RT));

NE_4 = exp( p.F * dPsi_m / p.RT);

J = x_C4 * OX_4 * (E_4 * cred_i - cox_i * NE_4);

end