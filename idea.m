clear
clc
close all

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

n_a = 200;
n_h = 20;
n_z = 11;
a_grid = linspace(0,40,n_a)';
z_grid = linspace(0.5,1.5,n_z)';
h_grid = linspace(0,1,n_h)';

% Move to GPU
a_grid_gpu = gpuArray(a_grid);
z_grid_gpu = gpuArray(z_grid);
h_grid_gpu = gpuArray(h_grid);

tic
ReturnMatrix = f_slow(a_grid_gpu,z_grid_gpu,h_grid_gpu,n_a,n_z);
toc 

tic
ReturnMatrix_fast = f_fast(a_grid_gpu,z_grid_gpu,h_grid_gpu,n_a,n_z);
toc 

err = max(abs(ReturnMatrix-ReturnMatrix_fast),[],"all")
