% First test case, solve the following equation
%                        ∂u/∂t - μ∆u = f,
% where f(t;μ) = 4μ(1-cos(t)), in 2D using P2 finite elements. 
% Exact solution is known to be 
%          u(x,t;μ) = 4 - ((x-0.5)^2 + (y-0.5).^2 + 4μ*sin(t))

clear all
close all
clc
% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;
t0 = 0; T = 1; dt = 0.1;
n2 = 25;
n1 = n2*L;

% Create and display the mesh
mesh = create_mesh(0,0,L,H,n1,n2);
%draw_mesh(mesh);

% Set parameter μ, function f and boundary conditions
mu = @(x) 5;
fun = @(t,x) 4*mu(x)*(1-cos(t))*x(1,:).^0;
u0 = @(x) 4 - ((x(1)-0.5).^2 + (x(2)-0.5).^2);
d1 = @(t,x) 4 - ((x(1)-0.5).^2 + (-0.5).^2 + 4*mu(x)*sin(t));
d2 = @(t,x) 4 - ((1-0.5).^2 + (x(2)-0.5).^2 + 4*mu(x)*sin(t));
d3 = @(t,x) 4 - ((x(1)-0.5).^2 + (1-0.5).^2 + 4*mu(x)*sin(t));
d4 = @(t,x) 4 - ((0-0.5).^2 + (x(2)-0.5).^2 + 4*mu(x)*sin(t));

% Derivative of the Dirichlet BC w.r.t time, used to apply lifting
d1_dt = @(t,x) -4*mu(x)*cos(t);
d2_dt = @(t,x) -4*mu(x)*cos(t);
d3_dt = @(t,x) -4*mu(x)*cos(t);
d4_dt = @(t,x) -4*mu(x)*cos(t);

dirichlet_functions = @(t,x) [d1(t,x); d2(t,x); d3(t,x); d4(t,x)]';
dirichlet_dt = @(t,x) [d1_dt(t,x); d2_dt(t,x); d3_dt(t,x); d4_dt(t,x)]';

% Create finite element space
bc = [1 1 1 1];

poly_degree = 'P2';
fespace = create_fespace(mesh,poly_degree,bc);

% options
opts.integrate_f = 1;

method = 'rk2_mid';

sol = solver_heat_rk(fespace,t0,T,dt,fun,u0,mu,dirichlet_functions, dirichlet_dt, method, opts);
% sol = solver_heat_BE(fespace,t0,T,dt,fun,u0,mu,dirichlet_functions, opts);

opts.namefile = 'data';
visualize_heat_solution(sol,0.1,opts)

%% Visualize exact solution
dt = 0.1;
t = linspace(t0,T, (T-t0)/dt + 1);
u_ex = @(t,x) 4 - ((x(1)-0.5).^2 + (x(2)-0.5).^2 + 4*mu(x)*sin(t));

u_ex_vec = zeros(length(fespace.nodes), length(t));

for j = 1:length(t)
    for i=1:length(fespace.nodes)
        u_ex_vec(i, j) = u_ex(t(j), sol.fespace.nodes(i, 1:2));
    end
end

sol_ex.u = u_ex_vec;
sol_ex.t0 = t0;
sol_ex.T = T;
sol_ex.dt = dt;
sol_ex.mesh = fespace.mesh;
sol_ex.fespace = fespace;

visualize_heat_solution(sol_ex,0.1,opts)

%% Convergence of the scheme
n = 5;
scale = 1;
dt_init = scale*0.1;
dts = [];
method = 'rk2_mid';

total_err_L2 = zeros(n,1);
total_err_H1 = zeros(n,1);

mu = @(x) 5;

% Set expression for exact solution
u_ex = @(t,x) 4 - ((x(1)-0.5).^2 + (x(2)-0.5).^2 + 4*mu(x)*sin(t));

% Set derivatives of the exact solution, used to compute error
u_dx = @(t,x) -2*(x(1)-0.5);
u_dy = @(t,x) -2*(x(2)-0.5);
u_dt = @(t,x) -4*mu(x)*cos(t);


for i = 1:n
    
    dt = dt_init/2^(i-1);
    dts = [dts;dt];
    
    sol = solver_heat_rk(fespace,t0,scale*T,dt,fun,u0,mu,dirichlet_functions, dirichlet_dt, method, opts);
%     sol = solver_heat_BE(fespace,t0,scale*T,dt,fun,u0,mu,dirichlet_functions, opts);

    % check error
    n_timesteps = size(sol.u,2);
    
    err_L2 = [];
    err_H1 = [];
    
    count = 0;
    t = t0;
    while(count < n_timesteps)
        count = count + 1;            
        
        exact_sol = u_ex(t, sol.fespace.nodes(1:length(fespace.nodes), 1:2));
        
        % Compute L2 and H1 error on the mesh for each timestep. Note that
        % compute_L2_error return the square of the error
        err_L2 = [err_L2, compute_L2_error(fespace, sol.u(:,count), @(x) u_ex(t,x))];
        err_H1 = [err_H1, compute_H1_error(fespace, sol.u(:,count),@(x) u_ex(t,x),@(x) [u_dx(t,x);u_dy(t,x)])];
        
        t = t + dt;
    end
    
    total_err_L2(i) = sqrt(sum(err_L2)*dt);
    total_err_H1(i) = sqrt(sum(err_H1)*dt);  
end

loglog(dts,total_err_L2,'.-r','Linewidth',1,'Markersize',20)
hold on

loglog(dts,total_err_H1,'.-b','Linewidth',1,'Markersize',20)

loglog(dts,1e+1*dts.^2,'--k')

title('L^2 and H^1 error');
legend('Global L^2 error on u','Global H^1 error on u','\tau^2');
