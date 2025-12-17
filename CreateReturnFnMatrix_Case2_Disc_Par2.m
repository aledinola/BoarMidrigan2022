function Fmatrix=CreateReturnFnMatrix_Case2_Disc_Par2(ReturnFn,n_da, n_a, n_z, d_gridvals, a_gridvals, z_gridvals, ReturnFnParams)
%This will shadow CreateReturnFnMatrix_Case2_Disc_Par2

% ALL: alpha,delta,r,tau,xi,iota,tau_a,xi_a,         tau_s,theta,gamma
% GROUP 1: alpha,delta,r,tau,xi,iota,tau_a,xi_a
% GROUP 2: tau_s,theta,gamma

n_d = n_da(1);

n_par1 = 8; % num. of parameters for f_ReturnFn1
ParamCell1=cell(n_par1,1);
for ii=1:n_par1
    ParamCell1(ii,1)={ReturnFnParams(ii)};
end

n_par2 = 3; % num. of parameters for f_ReturnFn2
ParamCell2=cell(n_par2,1);
for ii=1:n_par2 
    ParamCell2(ii,1)={ReturnFnParams(n_par1+ii)};
end

% ParamCell=cell(length(ReturnFnParams),1);
% for ii=1:length(ReturnFnParams)
%     if size(ReturnFnParams(ii))~=[1,1]
%         error('Using GPU for the return fn does not allow for any of ReturnFnParams to be anything but a scalar')
%     end
%     ParamCell(ii,1)={ReturnFnParams(ii)};
% end

if ~(isscalar(n_d) && isscalar(n_a) && isscalar(n_z))
    error('This function works only if d,a,z are all uni-dimensional')
end

d_grid = d_gridvals(1:n_d,1);
aprime_grid = reshape(a_gridvals,[1,n_a]);
a_grid      = reshape(a_gridvals,[1,1,n_a]);
z_grid      = reshape(z_gridvals,[1,1,1,n_z]);

%CashMatrix(h,1,a,z)
CashMatrix = arrayfun(@f_ReturnFn1,d_grid,a_grid,z_grid,ParamCell1{:});

%ReturnMatrix(h,a',a,z) = CashMatrix(h,1,a,z) + h_grid(h,1,1,1) + aprime_grid(1,a',1,1)
Fmatrix = arrayfun(@f_ReturnFn2,d_grid,aprime_grid,CashMatrix,ParamCell2{:});

end %end function


