% First test case, solve the following equation
%                        ∂u/∂t - ∆u = f,
% where f(x,t;μ) = 0.5*exp(-||x-xp||^2/(2sigma^2))*(1-cos(pi*t)), in 2D
% using P2 finite elements and Reduced Basis. 
% Exact solution is known to be 
%          u(x,t;μ) = 4 - ((x-0.5)^2 + (y-0.5).^2 + 4μ*sin(t))


clear all
close all
clc
% Set dimension of the domain and parameters of the mesh
L = 2;
H = 2;
t0 = 0; T = 0.5; dt = 0.05;
n2 = 25;
n1 = n2;

% Create and display the mesh
mesh = create_mesh(0,0,L,H,n1,n2);
%draw_mesh(mesh);

mu = @(x) 1;
xp = [1., 0.75];
sigma = 0.75;
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


n4 = @(t,x) 0;
%n4 = @(t,x) 2*(x(2)-0.5).*exp(-(-0.5).^2 - (x(2)-0.5).^2 -4*mu(x)*t);
dirichlet_functions = @(t,x) [d1(t,x); d2(t,x); d3(t,x); d4(t,x)]';
dirichlet_dt = @(t,x) [d1_dt(t,x); d2_dt(t,x); d3_dt(t,x); d4_dt(t,x)]';

neumann_functions = @(t,x) [0 ; 0; 0; n4(t,x)]';
neumann_dt = @(t,x) [0 ; 0; 0; n4(t,x)]';

% Create finite element space
bc = [1 1 1 1];

poly_degree = 'P2';
fespace = create_fespace(mesh,poly_degree,bc);

% options
opts.integrate_f = 1;
opts.integrate_neumann = 1;

method = 'rk2_mid';

sol = solver_heat_rk(fespace,t0,T,dt,fun,u0,mu,dirichlet_functions, dirichlet_dt, neumann_functions, neumann_dt, method, opts);
% sol = solver_heat_BE(fespace,t0,T,dt,fun,u0,mu,dirichlet_functions,neumann_functions, opts);

opts.namefile = 'data';
visualize_heat_solution(sol,0.5,opts)

%%

load('heat/sol_fine_2_scale1_L=2_4')


%%
n = 5;
scale = 1;
dt_init = scale*0.05;
dts = [];
method = 'rk2_mid';
total_err_L2 = zeros(n,1);
total_err_H1 = zeros(n,1);

mu = @(x) 1;

for i = 1:n
    
    dt = dt_init/2^(i-1);
    dts = [dts;dt];
    
    sol = solver_heat_rk(fespace,t0,scale*T,dt,fun,u0,mu,dirichlet_functions, dirichlet_dt, neumann_functions, neumann_dt, method, opts);
%     sol = solver_heat_BE(fespace,t0,scale*T,dt,fun,u0,mu,dirichlet_functions,neumann_functions, opts);

    % check error
    n_timesteps = size(sol.u,2);
    
%     bc_flags_fine = sol_fine.fespace.nodes(:,3:3);
%     bc_flags_coarse = sol.fespace.nodes(:,3:3);
    
%     idx_nodes_int_fine = find(bc_flags_fine==0);
%     idx_nodes_int_coarse = find(bc_flags_coarse==0);
    
    nodes_fine = sol_fine.fespace.nodes(:,1:2);
    nodes_coarse = sol.fespace.nodes(:,1:2);
%     nodes_fine = sol_fine.fespace.nodes(idx_nodes_int_fine,1:2);
%     nodes_coarse = sol.fespace.nodes(idx_nodes_int_coarse,1:2);
    t_fine = linspace(t0,scale*T, (scale*T-t0)/(dt_init/2^7) + 1);
    t_coarse = linspace(t0,scale*T, (scale*T-t0)/dt + 1);

    [~,idx_x] = ismembertol(nodes_coarse,nodes_fine,'ByRows', 1e-7);
    [~,idx_t] = ismembertol(t_coarse,t_fine, 1e-7);
    
    u_ref = sol_fine.u(idx_x,idx_t');

    
    err_L2 = [];
    err_H1 = [];
    
    count = 0;
    t = t0;
    while(count < n_timesteps)
        count = count + 1;    
        
        err_L2 = [err_L2, compute_L2_error_ref(fespace, sol.u(:,count), u_ref(:,count))];
        err_H1 = [err_H1, compute_H1_error_ref(fespace, sol.u(:,count), u_ref(:,count))];
        
        t = t + dt; 
    end
    
    total_err_L2(i) = sqrt(sum(err_L2)*dt);
    total_err_H1(i) = sqrt(sum(err_H1)*dt);  
end

loglog(dts,total_err_L2,'.-r','Linewidth',1,'Markersize',20)
hold on

loglog(dts,total_err_H1,'.-b','Linewidth',1,'Markersize',20)

loglog(dts,8e-2*dts.^2,'--k')

legend('Global L_2 error on u','Global H^1 error on u','\tau^2');

