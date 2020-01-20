% Second test case, solve the following equation
%                        ∂u/∂t - ∆u = f,
% where f(x,t;μ) = 0.5*exp(-||x-xp||^2/(2sigma^2))*(1-cos(pi*t)), in 2D
% using P2 finite elements and Reduced Basis. Exact solution is not known. 
% Convergence is studied using a reference solution.

clear all
close all
clc
% Set dimension of the domain and parameters of the mesh
L = 2;
H = 2;
scale = 1;
t0 = 0; T = scale*1/2; dt = scale*0.05;
n2 = 25;
n1 = n2;

% Create and display the mesh
mesh = create_mesh(0,0,L,H,n1,n2);
%draw_mesh(mesh);

% Create finite element space
bc = [1 1 1 1];

poly_degree = 'P2';
fespace = create_fespace(mesh,poly_degree,bc);

% options
opts.integrate_f = 1;

%% Load workspace with computed snapshots and RB. If executed, no need to run the three following cells

load('heat/data/workspace_RB_test2.mat')

%% Compute snapshots 


mu = @(x) 1;
xps = linspace(0.5,1.5,12);
sigmas = linspace(0.01,1.2,27);
sols = [];

for sigma = sigmas
    for i = 1:length(xps)
        for j = 1:length(xps)
            xp = [xps(i), xps(j)];
            
            fun = @(t,x) exp(-((x(1,:)-xp(1)).^2 + (x(2,:)-xp(2)).^2) / (2*sigma^2))*0.5*(1-cos(pi*t));
            u0 = @(x) 0*x(1);
            d1 = @(t,x) 0;
            d2 = @(t,x) 0;
            d3 = @(t,x) 0;
            d4 = @(t,x) 0;
            
            d1_dt = @(t,x) 0*x(1);
            d2_dt = @(t,x) 0*x(1);
            d3_dt = @(t,x) 0*x(1);
            d4_dt = @(t,x) 0*x(1);

            dirichlet_functions = @(t,x) [d1(t,x); d2(t,x); d3(t,x); d4(t,x)]';
            dirichlet_dt = @(t,x) [d1_dt(t,x); d2_dt(t,x); d3_dt(t,x); d4_dt(t,x)]';
    
            method = 'rk2_mid';
    
            sol = solver_heat_rk(fespace,t0,T,dt,fun,u0,mu,dirichlet_functions, dirichlet_dt, method, opts);
    
            % Take the interior dofs
            for k = 1:size(sol.u, 2)
                sol.u(:,k) = apply_dirichlet_bc_rhs(sol.u(:,k), fespace, @(x) [0;0;0;0]);
            end
            sols = [sols, sol];
        end
    end
end

%% Build matrix of snapshots

N_sol = length(sols);
S = [];
for i = 1:N_sol
    for t = 1: (T-t0)/dt + 1
        S = horzcat(S, sols(i).u(:, t));
    end
end


%% Build the reduced basis
n_nodes = length(fespace.nodes);

dt = scale*0.05;
S1 = S; 

V1 = POD(S1, dt);
disp('V1 done');

Proj = V1*V1';
S2 = (eye(n_nodes, n_nodes)-Proj)*S1;

V2 = POD(S2, dt);
disp('V2 done');

V = [V1, V2];

% M = assemble_mass(fespace);
% A = assemble_stiffness(mu, fespace);
% Xh = M + A;
% 
% V1 = POD_Xh(S1,Xh,dt^2);
% disp('V1 done');
% 
% Proj = V1*V1'*Xh;
% S2 = (eye(n_nodes, n_nodes)-Proj)*S1;
% 
% V2 = POD_Xh(S2, Xh, dt);
% disp('V2 done');
% 
% V = [V1, V2];


%% Test for an arbitrary parameter vector μ=(xp,sigma)
dt = 0.05;
mu = @(x) 1;

% Case b (refer to report, Fig.3)
xp = [1.25, 1.25];
sigma = 0.5;

%Case d (refer to report, Fig.3)
% xp = [1., 0.75];
% sigma = 0.75;

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

method = 'rk2_mid';

sol_RB = solver_heat_RB_rk(fespace,V,t0,T,dt,fun,u0,mu,dirichlet_functions, dirichlet_dt, method,opts);
% sol_RB = solver_heat_RB_rk_Xh(fespace,V,Xh,t0,T,dt,fun,u0,mu,dirichlet_functions, dirichlet_dt, method,opts);

visualize_heat_solution(sol_RB,0.3,opts)


%% Fine solution for % xp = [1.25, 1.25], sigma = 0.5

load('heat/data/sol_fine_test2_b_scale1_L=2')


%% Fine solution for xp = [1., 0.75], sigma = 0.75

load('heat/data/sol_fine_test2_d_scale1_L=2')

%% Convergence test with reference solution
n = 5;
scale = 1;
T = scale*1/2;
dt_init = scale*0.05;
dts = [];
method = 'rk2_mid';
total_err_L2 = zeros(n,1);
total_err_H1 = zeros(n,1);

mu = @(x) 1;

for i = 1:n
    
    dt = dt_init/2^(i-1);
    dts = [dts;dt];
    
%     sol_RB = solver_heat_RB_rk_Xh(fespace,V,Xh,t0,T,dt,fun,u0,mu,dirichlet_functions, dirichlet_dt, method,opts);
    sol_RB = solver_heat_RB_rk(fespace,V,t0,T,dt,fun,u0,mu,dirichlet_functions, dirichlet_dt, method,opts);   
    
    % Extract the nodes of the fine and coarse meshes
    nodes_fine = sol_fine.fespace.nodes(:,1:2);
    nodes_coarse = sol.fespace.nodes(:,1:2);

    % Extract the nodes of the fine and coarse temporal grid
    t_fine = linspace(t0,scale*T, (scale*T-t0)/(dt_init/2^7) + 1);
    t_coarse = linspace(t0,scale*T, (scale*T-t0)/dt + 1);

    % Find the indices of the nodes where the fine solution is computed on
    % the same nodes as the coarse solution
    [~,idx_x] = ismembertol(nodes_coarse,nodes_fine,'ByRows', 1e-10);
    [~,idx_t] = ismembertol(t_coarse,t_fine, 1e-10);
    
    % Take the value of the fine solution on the subgrid of the coarse
    % solution
    u_ref = sol_fine.u(idx_x,idx_t');

    err_L2 = [];
    err_H1 = [];
    
    count = 0;
    t = t0;
    
    n_timesteps = size(sol_RB.u,2);

    while(count < n_timesteps)
        count = count + 1;    
        
        % Compute error w.r.t reference solution
        err_L2 = [err_L2, compute_L2_error_ref(fespace, sol_RB.u(:,count), u_ref(:,count))];
        err_H1 = [err_H1, compute_H1_error_ref(fespace, sol_RB.u(:,count), u_ref(:,count))];
        
        t = t + dt; 
    end
    
    total_err_L2(i) = sqrt(sum(err_L2)*dt);
    total_err_H1(i) = sqrt(sum(err_H1)*dt);  
end

loglog(dts,total_err_L2,'.-r','Linewidth',1,'Markersize',20)
hold on

loglog(dts,total_err_H1,'.-b','Linewidth',1,'Markersize',20)

loglog(dts,5e-2*dts.^2,'--k')
title('L^2 and H^1 error');

legend('Global L^2 error on u','Global H^1 error on u','\tau^2');
