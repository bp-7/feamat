function [solution] = solver_heat_rk(fespace,t0,T,dt,fun,u0,mu,dirichlet_functions, dirichlet_dt, method,varargin)
integrate_f = 1;
if (nargin >= 11)
    opts = varargin{1};
    integrate_f = opts.integrate_f;
end

n_nodes = size(fespace.nodes,1);

% Load RK coefficients
coeffs = rk_coeffs(method);
nstages = coeffs.nstages;

% Compute initial condition
u = project_function(fespace, @(x) u0(x)');

% assemble constant matrices
A = assemble_stiffness(mu,fespace);
M = assemble_mass(fespace);

% Apply Dirichlet boundary condition
A = apply_dirichlet_bc_matrix(A,fespace,0);
M = apply_dirichlet_bc_matrix(M,fespace,1);

t = t0;
sol = zeros(n_nodes,ceil((T-t0)/dt));
sol(:,1) = u;

count = 1;
minnorm = min(sqrt(u.^2));
maxnorm = max(sqrt(u.^2));
minp = inf;
maxp = 0;

% Get interior dofs
u_int = apply_dirichlet_bc_rhs(u,fespace,@(x) [0;0;0;0]);

stages = zeros(length(u),nstages);
    
while (T-t>dt/2)
    count = count + 1;
    
    disp(['Time = ',num2str(t+dt)]);

    
    for j = 1:nstages

        cj = sum(coeffs.alpha(j,:));
                
        % Handle rhs
        fun_f = @(x) fun(t+cj*dt,x);
        dir = @(x) dirichlet_functions(t+cj*dt,x)';
        dir_dt = @(x) dirichlet_dt(t+cj*dt,x)';
        
        % Compute source term 
        f = zeros(n_nodes,1);
        ug_dt = zeros(n_nodes,1);
        ug = zeros(n_nodes,1);
   
        if (integrate_f)
            f = assemble_rhs(fespace,fun_f);
        end
        
        f = apply_dirichlet_bc_rhs(f,fespace,@(x) [0;0;0;0]);
        ug_dt = apply_dirichlet_bc_rhs(ug_dt,fespace,dir_dt);
        ug = apply_dirichlet_bc_rhs(ug,fespace,dir);
        
        if j == 1
            rhs = f - A*u_int - A*ug - M*ug_dt;
        else
            rhs = f - A*(u_int + ug + dt*(coeffs.alpha(j,1:j-1)*stages(:,1:j-1)')') - M*ug_dt;
        end
        rhs = apply_dirichlet_bc_rhs(rhs,fespace,@(x) [0;0;0;0]);
        
        % Compute the reduced LHS matrix
        mat = sparse(M + dt*coeffs.alpha(j,j)*A);
        
        stages(:,j) = mat\rhs;
    end
    
    % Compute boundary nodes at t=t_{n+1}
    dir = @(x) dirichlet_functions(t+dt,x)';
    ug = zeros(n_nodes,1);
    ug = apply_dirichlet_bc_rhs(ug,fespace,dir);
    
    % Apply RK update for interior dofs
    u_int = u_int + dt * stages * coeffs.b;
    
    % Add the boundary dofs to the interior dofs to recover u_{n+1}
    u = u_int + ug;  

    sol(:,count) = u;
    
    t = t + dt;
    
    minnorm = min(min(sqrt(u.^2)),minnorm);
    maxnorm = max(max(sqrt(u.^2)),maxnorm);
    minp = min(min(u),minp);
    maxp = max(max(u),maxp);
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
