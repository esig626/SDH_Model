function J = nhe1c(gNHE1,h_c,p)

% Taken from Falkenberg et al. (2010)

na_c = 15;
na_e = 135;
Kna = 0.00005;
Kh = 0.00005;
h_e = 0.00003981;
h_c = h_c * 10^(3);
J = gNHE1 * (na_e/(na_e + Kna) .* (h_c./(h_c + Kh)).^2 ...
                        - na_c./(na_c + Kna) .* (h_e./(h_e + Kh)).^2 );

end