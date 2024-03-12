function  J = akgd(gAKGD,akg_m,rsn_m,psn_m,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AKGDH + SUCLA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = gAKGD * ((akg_m ./ ( akg_m + p.KaKG )) ...
              .* (rsn_m ./( p.KRSn_aKG + rsn_m )) ...
              .* (psn_m./(p.KPSn_SCA + psn_m)));
          
          %  * p.coa / (p.Kcoa + p.coa)
end