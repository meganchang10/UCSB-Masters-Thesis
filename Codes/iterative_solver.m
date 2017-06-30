function u = iterative_solver(u,rhs,H,phi,dx,dy,m,n,iterations)

% This iterative solver is used in conjunction with our extrapolation
% methods, both third order and constant extrapolation. For third order, we
% solve 3 PDEs, and for first order we use only 1 PDE. The iterative solver
% uses a fictitious time step dtau to reach a steady state solution

% We begin by finding the first order one-sided first derivatives of our
% value u (instead of using WENO which is higher order and thus more 
% accurate, but this added accuracy is unnecessary in an iterative process,
% and using the WENO schemes is costly in time efficiency. Thus, we prefer
% the first order first derivatives which we can compute on the fly, and
% thus reduce memory costs as well.

% Important NOTE: The time step here is a fictitious time step. Thus it
% does not have the same time step restriction as our Stefan solver which
% is dependent on the advection dt < 0.5*dx/u. Instead the time step
% restriction we face is dt = 0.5 min(dx,dy) which will generally be much
% larger than the advection time step and allow us to extrapolate in fewer
% iterations

% Also, note, if the Heaviside value is equal to 0, the iterative solver
% will not change the value because a 0 H value means the value is known at
% that point. 

u_np1 = zeros(m,n);
u_np2 = zeros(m,n);
[nx,ny] = N_vector_generator(phi,m,n,dx,dy);
dtau = 0.5*min(dx,dy);

for t = 1:iterations
    for j = 2:n-1
        for i = 2:m-1
            ux_p = (u(i+1,j)-u(i,j))/dx;
            ux_m = (u(i,j)-u(i-1,j))/dx;
            uy_p = (u(i,j+1)-u(i,j))/dy;
            uy_m = (u(i,j)-u(i,j-1))/dy;
            
            if H(i,j) == 1
                nx_p = max(nx(i,j),0);
                nx_m = min(nx(i,j),0);
                ny_p = max(ny(i,j),0);
                ny_m = min(ny(i,j),0);
                u_np1(i,j) = u(i,j) + dtau*rhs(i,j) - dtau*H(i,j)*(nx_p*ux_m + nx_m*ux_p + ny_p*uy_m + ny_m*uy_p);
            elseif H(i,j) == 0
                u_np1(i,j) = u(i,j);
            end
        end
    end
    
    for j = 2:n-1
        for i = 2:m-1
            ux_p = (u_np1(i+1,j)-u_np1(i,j))/dx;
            ux_m = (u_np1(i,j)-u_np1(i-1,j))/dx;
            uy_p = (u_np1(i,j+1)-u_np1(i,j))/dy;
            uy_m = (u_np1(i,j)-u_np1(i,j-1))/dy;
            
            if H(i,j) == 1
                nx_p = max(nx(i,j),0);
                nx_m = min(nx(i,j),0);
                ny_p = max(ny(i,j),0);
                ny_m = min(ny(i,j),0);
                u_np2(i,j) = u_np1(i,j) + dtau*rhs(i,j) - dtau*H(i,j)*(nx_p*ux_m + nx_m*ux_p + ny_p*uy_m + ny_m*uy_p);
            elseif H(i,j) == 0
                u_np2(i,j) = u(i,j);
            end
        end
    end
    u = (u + u_np2)/2;
end

end