function J = KH(x_KH,k_m,h_m,p)

J = x_KH * (k_m * p.H_i - p.K_i * h_m );

end