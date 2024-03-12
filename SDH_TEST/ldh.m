function J = ldh(gLDH,pyr_c,nadh_c,nad_c,lac_c,p,h_c)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lactate Dehydrogenase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kpyr = p.Kpyr/2;
KH = 0.00005;
h_cc = h_c *10^3;
J = gLDH .* ( (pyr_c./(pyr_c + Kpyr)).^2 ...
    .* (nadh_c./( nadh_c + p.KRSp_PYR)) ...
    .* (h_cc./( h_cc + KH)).^2 ...
    - (lac_c./(lac_c + Kpyr*1000)).^2  ...
    .*   (nad_c./( nad_c + p.KRSp_PYR)) );

end