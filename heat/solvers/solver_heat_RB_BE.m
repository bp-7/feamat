function [solution] = solver_heat_RB_BE(fespace,V,t0,T,dt,fun,u0,mu,dirichlet_functions,varargin)
integrate_f = 1;
if (nargin >= 10)
    opts = varargin{1};
    integrate_f = opts.integrate_f;
end

n_nodes = size(fespace.nodes,1);

% Set initial condition
u = project_function(fespace, @(x) u0(x)');

% Assemble constant matrices
A = assemble_stiffness(mu,fespace);
M = assemble_mass(fespace);

% A = apply_dirichlet_bc_matrix(A,fespace,1);
M = apply_dirichlet_bc_matrix(M,fespace,1);

% Computation of reduced matrices
M_N = V'*M*V;
A_N = V'*A*V;

t = t0;
sol = zeros(n_nodes,ceil((T-t0)/dt));
sol(:,1) = u;

count = 1;
minnorm = min(sqrt(u.^2));
maxnorm = max(sqrt(u.^2));
minp = inf;
maxp = 0;

% Get interior dofs and project onto V
u_int = apply_dirichlet_bc_rhs(u, fespace, @(x) [0;0;0;0]);
u_int_N = V'*u_int;

while (T-t>dt/2)
    count = count + 1;
    t = t + dt;
    
    disp(['Time = ',num2str(t)]);
    
    fun_f = @(x) fun(t,x);
    dir = @(x) dirichlet_functions(t+dt,x)';
    
    if (integrate_f)
        f = assemble_rhs(fespace,fun_f);
    end
    
    % Compute reduced vector function
    f_N = V'*f;
    
    % Compute reduced LHS matrix
    mat_N = 1/dt * M_N + A_N;
    
    % Compute reduced RHS using lifting method
    rhs_N = f_N + 1/dt * M_N * u_int_N;
  
    % Solve for the interior dofs on the reduced basis and project onto high-fidelity space 
    u_int_N = mat_N\rhs_N;
    
    u_int = V*u_int_N;
    
    % Compute boundary dofs at t=t_{n+1}
    ug = zeros(n_nodes, 1);
    ug = apply_dirichlet_bc_rhs(ug,fespace,dir);
    
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
