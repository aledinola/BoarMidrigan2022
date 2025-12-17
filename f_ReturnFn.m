function F=f_ReturnFn(h,aprime,a,z,alpha,delta,r,tau,xi,iota,tau_a,xi_a,tau_s,theta,gamma)
% Action space: (h,aprime,a,z), where
% INPUT   TOOLKIT                         DESCRIPTION
% h:      decision variable               Labor supply
% aprime: next-period endogenous state    Future assets
% a:      endogenous state                Current assets
% z:      exogenous stochastic state      Idiosyncratic labor ability

w = f_prices(r,alpha,delta);

pretax_income = w*h*z+r*a;

income_tax = pretax_income-((1-tau)/(1-xi))*(pretax_income^(1-xi)) - iota;
wealth_tax = a-((1-tau_a)/(1-xi_a))*(a^(1-xi_a));

% Cash on hand
resources = (pretax_income-income_tax) + (a-wealth_tax) - aprime;

% consumption
c = resources/(1+tau_s); 

% Infeasible consumption
if c <= 0
    F = -Inf;
    return
end

if theta==1
    utility_c = log(c);
else
    utility_c = c^(1-theta)/(1-theta);
end

disutility_h = (h^(1+gamma))/(1+gamma);
F = utility_c - disutility_h;

end %end function