function  J = pdh(gPDH,pyr_m,rsn_m,p)

J = gPDH .* ( pyr_m .* p.coa_m ./(pyr_m .* p.coa_m + p.Kpyrm))...
    * (rsn_m./(rsn_m + p.KRSn_PYRCoA)) ;
end