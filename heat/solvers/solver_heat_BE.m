function [solution] = solver_heat_BE(fespace,t0,T,dt,fun_f,u0,mu,dirichlet_functions,varargin)

integrate_f = 1;
if (nargin >= 10)
    opts = varargin{1};
    integrate_f = opts.integrate_f;
end

n_nodes = size(fespace.nodes,1);

u = zeros(n_nodes, 1);

% Set initial condition
u = project_function(fespace, @(x) u0(x)');

% Assemble constant matrices
A = assemble_stiffness(mu,fespace);
M = assemble_mass(fespace);

A = apply_dirichlet_bc_matrix(A,fespace,0); 
M = apply_dirichlet_bc_matrix(M,fespace,dt);

t = t0;
sol = zeros(n_nodes,ceil((T-t0)/dt));
sol(:,1) = u;

count = 1;
minnorm = min(sqrt(u.^2));
maxnorm = max(sqrt(u.^2));
minp = inf;
maxp = 0;

while (T-t>dt/2)
    count = count + 1;
    t = t + dt;
    
    disp(['Time = ',num2str(t)]);
    
    fun = @(x) fun_f(t,x);
    dir = @(x) dirichlet_functions(t,x)';
    dir_b = @(x) dirichlet_functions(t-dt,x)';
    
    if (integrate_f)
        f = assemble_rhs(fespace,fun);
    end

    % Compute boundary dofs at t=t_n
    ug_b = zeros(n_nodes, 1);
    ug_b = apply_dirichlet_bc_rhs(ug_b,fespace,dir_b);
    
    % Compute boundary dofs at t=t_{n+1}
    ug = zeros(n_nodes, 1);
    ug = apply_dirichlet_bc_rhs(ug,fespace,dir);
    
    % Compute interior dofs
    u_int = apply_dirichlet_bc_rhs(u,fespace,@(x) [0;0;0;0]);
    
    % Compute LHS matrix
    mat = 1./dt * M + A;
    mat = apply_dirichlet_bc_matrix(mat,fespace,1);
    
    % Compute RHS using lifting method
    rhs = 1./dt * M*(u_int + ug_b - ug) - A*ug + f;
    rhs = apply_dirichlet_bc_rhs(rhs,fespace,@(x) [0;0;0;0]);
    
    % Solve for the interior dofs
    u_int = mat\rhs;
    
    % Add the boundary dofs to the interior dofs to recover u_{n+1}
    u = u_int + ug;
        
    
    minnorm = min(min(sqrt(u.^2)),minnorm);
    maxnorm = max(max(sqrt(u.^2)),maxnorm);
    minp = min(min(u),minp);
    maxp = max(max(u),maxp);
    
    sol(:,count) = u;
end

solution.u = sol;
solution.t0 = t0;
solution.T = T;
solution.dt = dt;
solution.mesh = fespace.mesh;
solution.fespace = fespace;
solution.minnorm = minnorm;
solution.maxnorm = maxnorm;
solution.minp = minp;
solution.maxp = maxp;
