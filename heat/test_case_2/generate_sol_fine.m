% This script generate solution on fine mesh and fine temporal grid for the 
% following problem,
%                        ∂u/∂t - ∆u = f,
% where f(x,t;μ) = 0.5*exp(-||x-xp||^2/(2sigma^2))*(1-cos(pi*t)), in 2D
% using P2 finite elements and given Runge-Kutta method. Save the solution. 

clear all
close all
clc

% Set dimension of the domain and parameters of the mesh
L = 2;
H = 2;
scale = 1;
t0 = 0; T = scale*1/2; dt = scale*0.05/2^7;
n2 = 100;
n1 = n2;

% Create and display the mesh
mesh = create_mesh(0,0,L,H,n1,n2);
%draw_mesh(mesh);

mu = @(x) 1;
xp = [1.25, 1.25];
sigma = 0.25;
fun = @(t,x) exp(-((x(1,:)-xp(1)).^2 + (x(2,:)-xp(2)).^2) / (2*sigma^2))*0.5*(1-cos(pi*t));
u0 = @(x) 0*x(1);
d1 = @(t,x) 0;
d2 = @(t,x) 0;
d3 = @(t,x) 0;
d4 = @(t,x) 0;

d1_dt = @(t,x) 0;
d2_dt = @(t,x) 0;
d3_dt = @(t,x) 0;
d4_dt = @(t,x) 0;


dirichlet_functions = @(t,x) [d1(t,x); d2(t,x); d3(t,x); d4(t,x)]';
dirichlet_dt = @(t,x) [d1_dt(t,x); d2_dt(t,x); d3_dt(t,x); d4_dt(t,x)]';

% Create finite element space
bc = [1 1 1 1];

poly_degree = 'P2';
fespace = create_fespace(mesh,poly_degree,bc);

% options
opts.integrate_f = 1;

method = 'rk2_mid';

sol_fine = solver_heat_rk(fespace,t0,T,dt,fun,u0,mu,dirichlet_functions, dirichlet_dt, method, opts);
% sol = solver_heat_BE(fespace,t0,T,dt,fun,u0,mu,dirichlet_functions, opts);

sol_fine.xp = xp;
solf_fine.sigma = sigma;

%visualize_heat_solution(sol_fine,0.05,opts)

save('heat/data/sol_fine_scale1', 'sol_fine')