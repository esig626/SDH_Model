function J = mct4(gMCT4,lac_c,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCT4 defined positive from inside to outside. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


KLAC = 1200;
Mce_lac = p.Mce_lac;
J = gMCT4 * ((lac_c./(lac_c + KLAC)) - (p.lac_e / (p.lac_e + Mce_lac)));

end