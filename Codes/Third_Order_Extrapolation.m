function u = Third_Order_Extrapolation(phi,u,m,n,dx,dy,iterations)

% Third order extrapolation uses the Heaviside functions which keeps track
% of which points can be considered knowns (i.e. when phi is negative). We
% solve 3 PDEs, and impose the steady state solution we want with the rhs
% input. We then use the iterative solver to iterate until we come to a
% steady state solution for the derivative input. For instance in the first
% round, we are tyring to solve for a second derivative u_nn with a rhs
% equal to zeros which means we will have constant extrapolation for our
% value u_nn.

% Then, we can use this second order derivative to impose our solution for
% the first derivative u_n. And use the first derivative u_n to impose a
% solution for u.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Important things to note:

% 1) When theta is very small, that is to say when the interface is too
% close to a grid point (theta < dx), we reassign the phi at that grid
% point to the interface so that phi_{i,j} = 0. Thus, it is considered
% a known value in our diffusion program since it cannot be contained in
% the vector of unknowns. HOWEVER, the value u at the location is not truly
% a known. So, when solving for our Heaviside functions which determine if
% we KNOW a value, we do not include 0. We do not have H_0(i,j) = 0 for 
% phi(i,j) <= 0, we simply have it for phi(i,j) < 0. In this way, we allow
% for the third order extrapolation method to still extrapolate a value for
% this point which will correct the value there and produce more accurate
% results

% 2) Our time step, dt is independent of the time step used in the Stefan
% solver because it does not depend on the same time step 

% 3) More on the HEAVISIDE FUNCTIONS: Since we use central differencing to
% calculate second-order first derivatives, we only consider the value of a
% derivative known if we have values for the points neighboring it on
% either side. Thus, for the first Heaviside function H(u,phi), which we
% named H_0 we can consider  it known (H(u,phi) = 0) if phi < 0. But to
% consider H(u_n,phi) known, we need to have known values at each of the 4
% neighboring points: u(i-1,j), u(i+1,j), u(i,j-1), u(i,j+1). In other
% words H_0 must be known at all these points. And then finally for
% H(u_nn,phi), we are looking for the second derivative or the first
% derivative of the first derivative. That is to say if we have known first
% derivatives at each of the neighboring points:  u_n(i-1,j), u_n(i+1,j),
% u_n(i,j-1), u_n(i,j+1), then we can consider the second derivative known
% and thus H(u_nn,phi) = 0. We assign the Heaviside function to zero
% because in the iterative solver, if a grid point has a Heaviside value of
% 0, then the value at that point is maintained, and not iterated or
% changed.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_n   = zeros(m,n);
u_nn  = zeros(m,n);
H_0   = ones(m,n);
H_n   = ones(m,n);
H_nn  = ones(m,n);

[nx,ny] = N_vector_generator(phi,m,n,dx,dy);

for j = 2:n-1
    for i = 2:m-1
        if phi(i,j) <= 0
            u_x = (u(i+1,j  ) - u(i-1,j  ))/(2*dx);
            u_y = (u(i  ,j+1) - u(i  ,j-1))/(2*dy);
            u_n(i,j)  = nx(i,j)*u_x + ny(i,j)*u_y;
        end
    end
end

for j = 2:n-1
    for i = 2:m-1
        if phi(i,j) <= 0
            temp_x = (u_n(i+1,j  ) - u_n(i-1,j  ))/(2*dx);
            temp_y = (u_n(i  ,j+1) - u_n(i  ,j-1))/(2*dy);
            u_nn(i,j)  = ny(i,j)*temp_y + nx(i,j)*temp_x;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solve for H's
for j = 2:n-1
    for i = 2:m-1
        if phi(i,j) < 0
            H_0(i,j) = 0;
        end
    end
end

for j = 2:n-1
    for i = 2:m-1
        if H_0(i,j-1) == 0 && H_0(i,j+1) == 0 && H_0(i-1,j) == 0 && H_0(i+1,j) == 0
            H_n(i,j) = 0;
        end
    end
end

for j = 2:n-1
    for i = 2:m-1
        if H_n(i,j-1) == 0 && H_n(i,j+1) == 0 && H_n(i-1,j) == 0 && H_n(i+1,j) == 0
            H_nn(i,j) = 0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhs = zeros(m,n);
u_nn = iterative_solver(u_nn,rhs,H_nn,phi,dx,dy,m,n,iterations);
rhs  = H_n.*u_nn;
u_n  = iterative_solver(u_n ,rhs,H_n ,phi,dx,dy,m,n,iterations);
rhs  = H_0.*u_n;
u    = iterative_solver(u   ,rhs,H_0 ,phi,dx,dy,m,n,iterations);

end