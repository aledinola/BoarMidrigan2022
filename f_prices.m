function w = f_prices(r,alpha,delta)

K_to_L = (alpha/(r+delta))^(1/(1-alpha));
w      = (1-alpha)*(K_to_L^alpha);

end %end function 