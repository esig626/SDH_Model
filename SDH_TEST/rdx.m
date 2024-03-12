function J = rdx(gRDX,nadh_c, nad_c,p,h_c)
h_cc = h_c .* 10^(3);
Kh = 0.00005;
J = gRDX .* ((nadh_c./(nadh_c + p.KRSn_2GAP))  .*  (h_cc./(h_cc + Kh*0.03))  ...
                                           - (nad_c./( nad_c + p.KRSp_PYR)));
end