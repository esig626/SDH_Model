function J = icd(gICD, cit_m, rsn_m,p) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Isocitrate dehydrogenase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = gICD * (cit_m /(cit_m + p.KCit)) * (rsn_m/(p.KRSn_CIT + rsn_m));
end