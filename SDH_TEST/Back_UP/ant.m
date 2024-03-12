function J = ant(x_ANT,adp_c,atp_c,Psi_i,atp_m,adp_m,Psi_x,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adenine Translocase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J= x_ANT * (adp_c / (adp_c + atp_c * exp(- p.F * Psi_i / p.RT))...
                - atp_m / (atp_m + adp_m * exp(- p.F* Psi_x / p.RT) ) );

end