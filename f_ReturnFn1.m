function cash = f_ReturnFn1(h,a,z,alpha,delta,r,tau,xi,iota,tau_a,xi_a)

w = f_prices(r,alpha,delta);

pretax_income = w*h*z+r*a;
income_tax = pretax_income-((1-tau)/(1-xi))*(pretax_income^(1-xi)) - iota;
wealth_tax = a-((1-tau_a)/(1-xi_a))*(a^(1-xi_a));
% Cash on hand
cash = (pretax_income-income_tax) + (a-wealth_tax);

end %end f_cash