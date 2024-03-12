function J = mdh(gMDH,mal_m,rsn_m,rpn_c,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Malate Dehydrogenase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J =  gMDH * ( mal_m ./(p.KMal + mal_m)) .* (rsn_m ./ (p.KRSn_MAL + rsn_m))...
                                            .* rpn_c./ (p.KRSn_MAL + rpn_c);

end