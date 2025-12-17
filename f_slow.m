function ReturnMatrix = f_slow(a_grid,z_grid,h_grid,n_a,n_z)

alpha = 0.33;
delta = 0.06;
r = 0.03;
tau = 0.26;
xi = 0.05;
iota = 0.16;
tau_a = 0;
xi_a = 0;
tau_s = 0.05;
theta = 1;
gamma = 2;

%h_grid
aprime_grid = reshape(a_grid,[1,n_a]);
a_grid      = reshape(a_grid,[1,1,n_a]);
z_grid      = reshape(z_grid,[1,1,1,n_z]);

%ReturnMatrix = zeros(n_h,n_a,n_a,n_z);

ReturnMatrix = arrayfun(@f_ReturnFn,h_grid,aprime_grid,a_grid,z_grid,r,tau_s,tau,xi,tau_a,xi_a,iota,theta,gamma,delta,alpha);

% for z_c=1:n_z
%     for a_c=1:n_a
%         for ap_c=1:n_a
%             for h_c=1:n_h
% 
%                 z = z_grid(z_c);
%                 a = a_grid(a_c);
%                 aprime = a_grid(ap_c);
%                 h = h_grid(h_c);
% 
%                 w = f_prices(r,alpha,delta);
% 
%                 pretax_income = w*h*z+r*a;
% 
%                 income_tax = pretax_income-((1-tau)/(1-xi))*(pretax_income^(1-xi)) - iota;
%                 wealth_tax = a-((1-tau_a)/(1-xi_a))*(a^(1-xi_a));
% 
%                 % Cash on hand
%                 resources = (pretax_income-income_tax) + (a-wealth_tax) - aprime;
% 
%                 % consumption
%                 c = resources/(1+tau_s);
% 
%                 F = -Inf;
%                 if c>0
%                     if theta==1
%                         utility_c = log(c);
%                     else
%                         utility_c = c^(1-theta)/(1-theta);
%                     end
% 
%                     disutility_h = (h^(1+gamma))/(1+gamma);
%                     F = utility_c - disutility_h;
%                 end
%                 ReturnMatrix(h_c,ap_c,a_c,z_c) = F;
% 
%             end %end h
%         end %end a'
%     end %end a
% end %end z


end %end function