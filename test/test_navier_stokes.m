clear all
close all
clc
% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

n2 = 20;
n1 = 20;

% Create and display the mesh
mesh = create_mesh(0,0,L,H,n1,n2);

f = @(t,x) [1;0];
nu = @(x) 1;
dirichlet_functions = @(t,x) [0 0;0 0;0 0;0 0]';
neumann_functions = @(t,x) [0 0;0 0;0 0;cos(2*pi*t)*0 0]';

% Create finite element space
bc = [1 0 1 0];

poly_degree = 'P2';
fespace_u = create_fespace(mesh,poly_degree,bc);

poly_degree = 'P1';
fespace_p = create_fespace(mesh,poly_degree,bc);

% options
opts.integrate_f = 1;
opts.integrate_neumann = 1;

sol = solver_navier_stokes(fespace_u,fespace_p,0,1,0.05,f,@(x) [0;0], @(x) 0,nu,dirichlet_functions,neumann_functions,opts);

opts.print = 0;
opts.namefile = 'data';
% opts.print_only = 'U';
opts.print_only = 0;
visualize_stokes_solution(sol,0.01,opts)

%% check convergence in time

clear all
close all
clc
% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

n2 = 30;
n1 = 30;

% Create and display the mesh
mesh = create_mesh(0,0,L,H,n1,n2);

% Exact solution
uex = @(t,x) sin(t*pi)*(x(2)-L).^3;
vex = @(t,x) t*cos(x(1)*pi);

% derivatives of exact solution
udx = @(t,x) 0;
udy = @(t,x) 3*sin(t*pi)*(x(2)-L).^2;
vdx = @(t,x) -pi*t*sin(x(1)*pi);
vdy = @(t,x) 0;
udxdx = @(t,x) 0;
udydy = @(t,x) 6*sin(t*pi)*(x(2)-L);
vdxdx = @(t,x) -pi^2*t*cos(x(1)*pi);
vdydy = @(t,x) 0;

udt = @(t,x) pi*cos(t*pi)*(x(2)-L).^3;
vdt = @(t,x) cos(x(1)*pi);

w = 2;
pex = @(t,x) (L-x(1))/L*cos(t*w);
pdx = @(t,x) -cos(t*w)/L;
pdy = @(t,x) 0;

f = @(t,x) [udt(t,x) - udxdx(t,x) - udydy(t,x) + uex(t,x)*udx(t,x) + ...
    vex(t,x)*udy(t,x) + pdx(t,x);
    vdt(t,x) - vdxdx(t,x) - vdydy(t,x) + uex(t,x)*vdx(t,x) + ...
    vex(t,x)*vdy(t,x) + pdy(t,x)];
center_nu = [L/4;H/2]';
nu = @(x) 1;
dirichlet_functions = @(t,x) [uex(t,x) vex(t,x);
    uex(t,x) vex(t,x);
    uex(t,x) vex(t,x);
    uex(t,x) vex(t,x)]';
neumann_functions = @(t,x) [0 0;0 0;0 0;0 0]';

% Create finite element space
bc = [1 0 1 1];

poly_degree = 'P2';
fespace_u = create_fespace(mesh,poly_degree,bc);

poly_degree = 'P1';
fespace_p = create_fespace(mesh,poly_degree,bc);

% options
opts.integrate_f = 1;
opts.integrate_neumann = 1;


dt_init = 0.1;
t0 = 0;
T = dt_init*4;

n = 4;
dts = [];

total_err_ul2 = zeros(n,1);
total_err_uh1 = zeros(n,1);
total_err_pl2 = zeros(n,1);


for i = 1:n
    
    dt = dt_init/2^(i-1);
    dts = [dts;dt];
    
    sol = solver_navier_stokes(fespace_u,fespace_p,t0,T,dt,f,@(x) [uex(0,x);vex(0,x)],@(x) pex(0,x),nu,dirichlet_functions,neumann_functions,opts);
    
    % check error
    
    n_timesteps = size(sol.u,2);
    
    count = 0;
    errul2 = [];
    erruh1 = [];
    errpl2 = [];
    
    t = t0;
    while(count < n_timesteps)
        count = count + 1;
        err1 = compute_error(fespace_u,sol.u1(:,count),@(x) uex(t,x),@(x) [udx(t,x);udy(t,x)],'L2');
        err2 = compute_error(fespace_u,sol.u2(:,count),@(x) vex(t,x),@(x) [vdx(t,x);vdy(t,x)],'L2');
        errul2 = [errul2;err1^2+err2^2];
        err1 = compute_error(fespace_u,sol.u1(:,count),@(x) uex(t,x),@(x) [udx(t,x);udy(t,x)],'H1');
        err2 = compute_error(fespace_u,sol.u2(:,count),@(x) vex(t,x),@(x) [vdx(t,x);vdy(t,x)],'H1');
        erruh1 = [erruh1;err1^2+err2^2];
        err1 = compute_error(fespace_p,sol.p(:,count),@(x) pex(t,x),@(x) 0,'L2');
        errpl2 = [errpl2;err1^2];
        t = t + dt;
    end
        
    total_err_ul2(i) = sqrt(sum(errul2)*dt);
    total_err_uh1(i) = sqrt(sum(erruh1)*dt);
    total_err_pl2(i) = sqrt(sum(errpl2)*dt);
  
end

loglog(dts,total_err_ul2,'.-r','Linewidth',1,'Markersize',20)
hold on

loglog(dts,total_err_uh1,'.-b','Linewidth',1,'Markersize',20)
loglog(dts,total_err_pl2,'.-g','Linewidth',1,'Markersize',20)

loglog(dts,dts*1e-2,'--k')

legend('Global error on ul2','Global error on uh1','Global error on pl2','\Delta t');



%%
opts.print = 0;
opts.namefile = 'data';
% opts.print_only = 'U';
opts.print_only = 0;
visualize_stokes_solution(sol,0.01,opts)