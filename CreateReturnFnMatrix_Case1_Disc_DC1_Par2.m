function Fmatrix=CreateReturnFnMatrix_Case1_Disc_DC1_Par2(ReturnFn, n_d, n_z, d_gridvals, aprime_grid, a_grid, z_gridvals, ReturnFnParams, Level)
%If there is no d variable, just input n_d=0 and d_grid=0

n_a = length(a_grid);
n_aprime = size(aprime_grid,1);

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

d_grid = d_gridvals;

a_grid      = reshape(a_grid,[1,1,n_a]);
z_grid      = reshape(z_gridvals,[1,1,1,n_z]);

%CashMatrix(h,1,a,z)
CashMatrix = arrayfun(@f_ReturnFn1,d_grid,a_grid,z_grid,ParamCell1{:});

aprime_grid = reshape(aprime_grid,[1,n_aprime,n_a,n_z]);

%ReturnMatrix(h,a',a,z) = CashMatrix(h,1,a,z) + h_grid(h,1,1,1) + aprime_grid(1,a',1,1)
Fmatrix = arrayfun(@f_ReturnFn2,d_grid,aprime_grid,CashMatrix,ParamCell2{:});



end %end function


