function ReturnMatrix = f_fast(a_grid,z_grid,h_grid,n_a,n_z)

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

%CashMatrix(h,1,a,z)
CashMatrix = arrayfun(@f_ReturnFn1,h_grid,a_grid,z_grid,alpha,delta,r,tau,xi,iota,tau_a,xi_a);

%ReturnMatrix(h,a',a,z) = CashMatrix(h,1,a,z) + h_grid(h,1,1,1) + aprime_grid(1,a',1,1)
ReturnMatrix = arrayfun(@f_ReturnFn2,h_grid,aprime_grid,CashMatrix,tau_s,theta,gamma);

end %end function

