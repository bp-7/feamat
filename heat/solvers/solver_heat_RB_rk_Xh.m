function [solution] = solver_heat_RB_rk_Xh(fespace,V,Xh,t0,T,dt,fun,u0,mu,dirichlet_functions, dirichlet_dt, method,varargin)
integrate_f = 1;
if (nargin >= 13)
    opts = varargin{1};
    integrate_f = opts.integrate_f;
end

% Load RK coefficients
coeffs = rk_coeffs(method);
nstages = coeffs.nstages;

n_nodes = size(fespace.nodes,1);
N = size(V,2);

% Compute initial condition
u = project_function(fespace, @(x) u0(x)');

% assemble constant matrices
A = assemble_stiffness(mu,fespace);
M = assemble_mass(fespace);

% Apply Dirichlet boundary condition
A = apply_dirichlet_bc_matrix(A,fespace,0);
M = apply_dirichlet_bc_matrix(M,fespace,1);

% Computation of reduced matrices
M_N = V'*Xh*M*V;
A_N = V'*Xh*A*V;

t = t0;
sol = zeros(n_nodes,ceil((T-t0)/dt));
sol(:,1) = u;

count = 1;
minnorm = min(sqrt(u.^2));
maxnorm = max(sqrt(u.^2));
minp = inf;
maxp = 0;

% Get interior dofs and project onto V
u_int = apply_dirichlet_bc_rhs(u,fespace,@(x) [0;0;0;0]);
u_int_N = V'*Xh*u_int;

% Parameter for the update of the solution, ={0,1}
option = 0;

while (T-t>dt)

    disp(['Time = ',num2str(t+dt)]);

    stages = zeros(N,nstages);
    
    u_int_N = V'*Xh*u_int;
        
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
        
        % Compute derivative of boundary dofs and boundary dofs at t=t_{n+1}
        ug_dt = apply_dirichlet_bc_rhs(ug_dt,fespace,dir_dt);
        ug = apply_dirichlet_bc_rhs(ug,fespace,dir);
        
        f_N = V'*Xh*f;
        ug_dt_N = V'*Xh*ug_dt;
           
        if j == 1
            rhs = f_N - A_N*u_int_N - V'*Xh*A*ug - ug_dt_N;
        else
            rhs = f_N - A_N*u_int_N - V'*Xh*A*ug - A_N*dt*coeffs.alpha(j,1:j-1).*stages(:,1:j-1) - ug_dt_N;
        end
        
        % Compute the reduced LHS matrix
        mat_N = M_N + dt*coeffs.alpha(j,j)*A_N;
       
        stages(:,j) = mat_N\rhs;
    end
    
    % Compute boundary nodes at t=t_{n+1}
    dir = @(x) dirichlet_functions(t+dt,x)';
    ug = zeros(n_nodes,1);
    ug = apply_dirichlet_bc_rhs(ug,fespace,dir);
        
    % Here two possibilities to update the solution on the interior dofs
    if option 
        u_int_N = u_int_N + dt * stages * coeffs.b;
        u_int = V*u_int_N;
    else
        u_int = u_int + dt * V * stages * coeffs.b;
    end
    
    % Add the boundary dofs to the interior dofs to recover u_{n+1}
    u = u_int + ug;
    
    count = count + 1;
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
