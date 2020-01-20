function visualize_heat_solution(solution, delay,varargin)
print = 0;
if (nargin >= 3)
    opts = varargin{1};
    namefile = opts.namefile;
end

u = solution.u;
t0 = solution.t0;
T = solution.T;
dt = solution.dt;
mesh = solution.mesh;
fespace = solution.fespace;
L = mesh.L;
H = mesh.H;
n_timesteps = size(u,2);
t = t0;
count = 0;
figure(1)

while (count < n_timesteps)  
    count = count + 1;
    plot_fe_function(u(:, count), fespace);
    title(['U at t = ',num2str(t)])
    axis([0 L 0 H])
    colorbar
    shading interp

    pbaspect([L H 1])
    pause(delay)
    if (print)
        saveas(gcf,['data/',namefile','_',num2str(count),'.png'])
    end
    t = t + dt;
end
    
end
        
       