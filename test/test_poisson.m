clear all
close all
clc

% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

n2 = 50;
n1 = n2*L;

% Create and display the mesh
mesh = create_mesh(L,H,n1,n2);
% draw_mesh(mesh);

solex = @(x,y) sin(pi*x).*sin(pi*y);

f = @(x) 2*pi^2*sin(pi*x(1))*sin(pi*x(2));
mu = @(x) 1;
dirichlet_functions = @(x) [x(1);0;0;0];
neumann_functions = @(x) [0;0;0;0];

% Create finite element space
bc = [1 0 0 0]; 

poly_degree = 'P1';
fespace = create_fespace(mesh,poly_degree,bc);

% Assemble matrix and rhs
[A,b] = assembler_poisson(f,mu,fespace,neumann_functions);

% Apply Dirichlet boundary conditions
[A,b] = apply_bc(A,b,fespace,dirichlet_functions);

% Solve the linear system
sol = A\b;

n1 = size(mesh.X,1);
n2 = size(mesh.X,2);

figure
surf(mesh.X,mesh.Y,reshape(sol,n1,n2),'EdgeColor','none','LineStyle','none','FaceLighting','phong');
pbaspect([1 1 1])


l2norm = compute_norm(fespace,sol,'L2');

display(['Norm = ', num2str(l2norm)]);

%% Here we check the convergence of the error
clear all
clc

% Set dimension of the domain and parameters of the mesh
L = 1;
H = 1;

solex = @(x) sin(pi*x(1)).*sin(pi*x(2));
gradex = @(x) [cos(pi*x(1)).*sin(pi*x(2));cos(pi*x(2)).*sin(pi*x(1))]*pi;


f = @(x) 2*pi^2*sin(pi*x(1))*sin(pi*x(2));
mu = @(x) 1;
dirichlet_functions = @(x) [0;0;0;0];

errl2 = [];
errh1 = [];
h = [];

for i = 1:5

    n2 = 5*2^(i-1);
    n1 = n2*L;

    % Create and display the mesh
    mesh = create_mesh(L,H,n1,n2);
    
    h = [h 1/n2];

    % Create finite element space
    bc = [1 1 1 1]; 

    poly_degree = 'P1';
    fespace = create_fespace(mesh,poly_degree,bc);

    % Assemble matrix and rhs
    [A,b] = assembler_poisson(f,mu,fespace);

    % Apply Dirichlet boundary conditions
    [A,b] = apply_bc(A,b,fespace,dirichlet_functions);

    % Solve the linear system
    sol = A\b;

    l2error = compute_error(fespace,sol,solex,gradex,'L2');
    l2norm = compute_norm(fespace,sol,'L2');

    disp(['Relative L2 error = ', num2str(l2error/l2norm)]);
    disp(' ');
    errl2 = [errl2 l2error];
    
    h1error = compute_error(fespace,sol,solex,gradex,'H1');
    h1norm = compute_norm(fespace,sol,'H1');

    disp(['Relative H1 error = ', num2str(h1error/h1norm)]);
    disp(' ');
    errh1 = [errh1 h1error];
    
end

loglog(h,errl2,'.-r','Linewidth',3,'Markersize',20)
hold on
loglog(h,h.^2,'--r','Linewidth',1);

loglog(h,errh1,'.-b','Linewidth',3,'Markersize',20)
hold on
loglog(h,h.^1,'--b','Linewidth',1);

minl2 = min(errl2);
minh1 = min(errh1);
maxl2 = max(errl2);
maxh1 = max(errh1);

legend('L2 error','h^2','H1 error','h','Location','Southeast');
pbaspect([1 1 1]);
axis([min(h) max(h) min([minl2 minh1]) max([maxl2 maxh1])])
set(gca,'Fontsize',25);

