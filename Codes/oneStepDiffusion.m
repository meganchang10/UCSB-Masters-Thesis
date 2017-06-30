function [u_n] = oneStepDiffusion(u_n,phi,m,n,D,S,BC_interface,BC_domain,BC_interface_old,x,y,dx,dy,dt,t)

% LEVEL SET: We use the level set method to keep track of our interface
% where the interface is the zero level set phi0, and the reacted region is
% phi < 0 and the unreacted region is phi > 0. 

% The diffusion solver uses a method called the Ghost-Fluid method which
% calculates the fractional distance theta from a point and the interface
% if the interface is contained between two neighboring points. This can be
% identified by checking for a sign change between two neighboring points
% phi such that phi(i,j)*phi(i+1,j) < 0 for instance. Since we know the 
% product would be negative if there was a sign change between the two
% points. In this case, we know the interface is contained between the two
% points and we calculate theta as the fractional distance of dx between
% phi(i,j) and the interface. Theta can then be used to calculate ghost
% points which are added to the RHS and account for the fraction of
% diffusion that will be felt from that side. However, theta is used in the
% denominator, so if it is too small, it will make the solution approach
% infinity and cause problems, so we reassign grid points with thetas that
% are too small to be on the interface, phi(i,j) = 0. We do this at the
% beginning of the diffusion code since in the next step we need to be able
% to count how many points are unknown, and points where phi = 0 are
% considered known by this solver since we can impose Dirichlet Boundary
% Conditions at the interface. 

for j = 2:n-1
    for i = 2:m-1
        if phi(i,j) == 0
            u_n(i,j) = BC_interface(i,j);
        elseif phi(i,j) < 0
            %LEFT

            if phi(i,j)*phi(i-1,j  ) < 0
                theta = abs(phi(i,j))/dx;
                if theta < dx
                    phi(i,j) = 0;
                    u_n(i,j) = BC_interface(i,j);
                    continue
                end
            end

            %RIGHT
            if phi(i,j)*phi(i+1,j  ) < 0
                theta = abs(phi(i,j))/dx;
                if theta < dx
                    phi(i,j) = 0;
                    u_n(i,j) = BC_interface(i,j);
                    continue
                end
            end
            
            %BELOW
            if phi(i,j)*phi(i,j-1) < 0
                theta = abs(phi(i,j))/dy;
                if theta < dy
                    phi(i,j) = 0;
                    u_n(i,j) = BC_interface(i,j);
                    continue
                end
            end
            
            %ABOVE
            if phi(i,j)*phi(i  ,j+1) < 0
                theta = abs(phi(i,j))/dy;
                if theta < dy
                    phi(i,j) = 0;
                    u_n(i,j) = BC_interface(i,j);
                    continue
                end
            end

        end
    end
end

% To build the most efficient solver, we build a sparse matrix since our
% diffusion matrix can only have at most 5 points per line, one for the
% diffusion contribution from the center and 4 neighboring points. This
% reduces the required memory storage necessary. The sparse matrices work
% using three vectors to track the value of the non-zero values A that 
% would be in our matrix, as well as their position i and j. I have named
% these vectors sparse_A, sparse_i, and sparse_j. To maintain efficiency,
% I need to preallocate their lengths by determining how many non-zero
% values we have. Thus, I filter through the entire domain once to count
% how many points have a neighboring point inside the reacted domain that
% will contribute diffusion. This can be found by seeing if the product of
% phi(i,j) and each of its neighboring points will be positive.

RHS_numel = 0;
A_numel   = 0;
nm = n*m;
tag = zeros(nm,1);
for j = 2:n-1
    for i = 2:m-1
        if phi(i,j) < 0
            RHS_numel = RHS_numel + 1;
            %Center
            A_numel = A_numel + 1;
            tag(p(i,j,m)) = RHS_numel;
            
            %LEFT
            if i > 2
                if phi(i,j)*phi(i-1,j  ) > 0
                    A_numel = A_numel + 1;
                end
            end
            
            %RIGHT
            if i < m-2
                if phi(i,j)*phi(i+1,j  ) > 0
                    A_numel = A_numel + 1;
                end
            end
            
            %BELOW
            if j > 2
                if phi(i,j)*phi(i  ,j-1) > 0
                    A_numel = A_numel + 1;
                end
            end
            
            %ABOVE 
            if j < n-2
                if phi(i,j)*phi(i  ,j+1) > 0
                    A_numel = A_numel + 1;
                end
            end
            
        end
    end
end

untag_i  = zeros(RHS_numel,1);
untag_j  = zeros(RHS_numel,1);
sparse_i = zeros(A_numel  ,1);
sparse_j = zeros(A_numel  ,1);
sparse_A = zeros(A_numel  ,1);
RHS      = zeros(RHS_numel,1);

% Here, we build our diffusion matrix A (in sparse form). When a point
% u(i,j) has a neighboring point that is also in the unreacted region, we
% calculate the diffusion coefficient and multiply it by 0.5 before storing
% it in sparse_A because we are implementing the implicit Crank-Nicholson 
% Scheme so that we can have a dt proportional to dx instead of dx^2. We
% also keep track of the location by adding i to sparese_i and j to 
% sparse_j. In the case where the interface is containted between the point
% i,j and the neighboring point, we calculate theta and use it in our
% formula for ghost values which will be added to the RHS of our equation
% later. For the Crank-Nicolson scheme, keep in mind that we solve using 
% values from the time t and t-dt. Thus, we have two ghost values for 
% ghost_nm1 and ghost_n. We do this for each neighboring point u(i-1,j),
% u(i+1,j), u(i,j-1), and u(i,j+1) 

count     = 0;
count_rhs = 0;

for j = 2:n-1
    for i = 2:m-1
        pp = p(i,j,m);
        tag_p = tag(pp);
        
        if phi(i,j) < 0
            %Reinitiate
            Below     = 0;
            Above     = 0;
            Right     = 0;
            Left      = 0;
            ghost_n   = 0;
            ghost_nm1 = 0;
            Center    = 2*D*dt/dy^2 + 2*D*dt/dx^2;
            
            %Diffusion from the grid point to the LEFT
            if i > 2
                if phi(i,j)*phi(i-1,j) > 0
                    Left  =                         -D*dt/dx^2;               %Left
                    count = count + 1;
                    sparse_i(count) = tag(p(i-1,j  ,m));
                    sparse_j(count) = tag_p;
                    sparse_A(count) = 0.5*Left; %Crank Nicolson Scheme
                elseif phi(i,j)*phi(i-1,j) < 0
                    theta = abs(phi(i,j))/dx;
                    interpolate_BC     = BC_interface(i,j)     + theta*(BC_interface(i-1,j)     -BC_interface(i,j));
                    interpolate_BC_old = BC_interface_old(i,j) + theta*(BC_interface_old(i-1,j) -BC_interface_old(i,j));
                    ghost_n   = ghost_n   + (1/theta)*(D*dt/dx^2)*interpolate_BC;
                    ghost_nm1 = ghost_nm1 + (1/theta)*(D*dt/dx^2)*interpolate_BC_old;
                    Center = Center + ((1-theta)/theta)*(D*dt/dx^2);
                elseif phi(i,j)*phi(i-1,j) == 0
                    ghost_n   = ghost_n   + D*dt/dx^2*BC_interface(i-1,j);
                    ghost_nm1 = ghost_nm1 + D*dt/dx^2*BC_interface_old(i-1,j);
                end
            else
                ghost_n   = ghost_n   + D*dt/dx^2*BC_domain;
                ghost_nm1 = ghost_nm1 + D*dt/dx^2*BC_domain;
            end

            %Diffusion from the grid point to the RIGHT
            if i < m-2
                if phi(i,j)*phi(i+1,j) > 0
                    Right  =                         -D*dt/dx^2;             %Right
                    count = count + 1;
                    sparse_i(count) = tag(p(i+1,j  ,m));
                    sparse_j(count) = tag_p;
                    sparse_A(count) = 0.5*Right;
                elseif phi(i,j)*phi(i+1,j) < 0
                    theta = abs(phi(i,j))/dx;
                    interpolate_BC     = BC_interface(i,j)     + theta*(BC_interface(i+1,j)     -BC_interface(i,j));
                    interpolate_BC_old = BC_interface_old(i,j) + theta*(BC_interface_old(i+1,j) -BC_interface_old(i,j));
                    ghost_n   = ghost_n   + (1/theta)*(D*dt/dx^2)*interpolate_BC;
                    ghost_nm1 = ghost_nm1 + (1/theta)*(D*dt/dx^2)*interpolate_BC_old;
                    Center = Center + ((1-theta)/theta)*(D*dt/dx^2);
                elseif phi(i,j)*phi(i+1,j) == 0
                    ghost_n   = ghost_n   + D*dt/dx^2*BC_interface(i+1,j);
                    ghost_nm1 = ghost_nm1 + D*dt/dx^2*BC_interface_old(i+1,j);
                end
            else
                ghost_n   = ghost_n   + D*dt/dx^2*BC_domain;
                ghost_nm1 = ghost_nm1 + D*dt/dx^2*BC_domain;
            end

            %Diffusion from the grid point BELOW
            if j > 2
                if phi(i,j)*phi(i,j-1) > 0
                    Below   =                         -D*dt/dy^2;             %Below
                    count = count + 1;
                    sparse_i(count) = tag_p;
                    sparse_j(count) = tag(p(i  ,j-1,m));
                    sparse_A(count) = 0.5*Below; %Crank Nicolson Scheme
                elseif phi(i,j)*phi(i,j-1) < 0
                    theta = abs(phi(i,j))/dy;
                    interpolate_BC     = BC_interface(i,j)     + theta*(BC_interface(i,j-1)     -BC_interface(i,j));
                    interpolate_BC_old = BC_interface_old(i,j) + theta*(BC_interface_old(i,j-1) -BC_interface_old(i,j));
                    ghost_n   = ghost_n   + (1/theta)*(D*dt/dy^2)*interpolate_BC;
                    ghost_nm1 = ghost_nm1 + (1/theta)*(D*dt/dx^2)*interpolate_BC_old;
                    Center = Center + ((1-theta)/theta)*(D*dt/dy^2);
                elseif phi(i,j)*phi(i,j-1) == 0
                    ghost_n   = ghost_n   + D*dt/dy^2*BC_interface(i,j-1);
                    ghost_nm1 = ghost_nm1 + D*dt/dy^2*BC_interface_old(i,j-1);
                end
            else
                ghost_n   = ghost_n   + D*dt/dy^2*BC_domain;
                ghost_nm1 = ghost_nm1 + D*dt/dy^2*BC_domain;
            end

            %Diffusion from the grid point ABOVE
            if j < n-2
                if phi(i,j)*phi(i,j+1) > 0
                    Above  =                         -D*dt/dy^2;             %Above
                    count = count + 1;
                    sparse_i(count) = tag_p;
                    sparse_j(count) = tag(p(i  ,j+1,m));
                    sparse_A(count) = 0.5*Above; %Crank Nicolson Scheme
                elseif phi(i,j)*phi(i,j+1) < 0
                    theta = abs(phi(i,j))/dy;
                    interpolate_BC     = BC_interface(i,j)     + theta*(BC_interface(i,j+1)    - BC_interface(i,j));
                    interpolate_BC_old = BC_interface_old(i,j) + theta*(BC_interface_old(i,j+1)- BC_interface_old(i,j));
                    ghost_n   = ghost_n   + (1/theta)*(D*dt/dy^2)*interpolate_BC;
                    ghost_nm1 = ghost_nm1 + (1/theta)*(D*dt/dx^2)*interpolate_BC_old;
                    Center = Center + ((1-theta)/theta)*(D*dt/dy^2);
                elseif phi(i,j)*phi(i,j+1) == 0
                    ghost_n   = ghost_n   + D*dt/dy^2*BC_interface(i,j+1);
                    ghost_nm1 = ghost_nm1 + D*dt/dy^2*BC_interface_old(i,j+1);
                end
            else
                ghost_n   = ghost_n   + D*dt/dy^2*BC_domain;
                ghost_nm1 = ghost_nm1 + D*dt/dy^2*BC_domain;
            end

            count     = count + 1;
            count_rhs = count_rhs + 1;
            
            untag_i(tag_p) = i;
            untag_j(tag_p) = j;
            sparse_i(count) = count_rhs;
            sparse_j(count) = count_rhs;
            sparse_A(count) = 0.5*Center + 1;
            
            %Crank Nicolsoning the u_n on the RHS side
            CL =   Left*u_n(i-1,j  );
            CR =  Right*u_n(i+1,j  );
            CB =  Below*u_n(i  ,j-1);
            CA =  Above*u_n(i  ,j+1);
            CC = Center*u_n(i  ,j  );
            Crank = CL + CR + CB + CA + CC;
            
            % Apply the Crank-Nicholson scheme which takes an average of
            % values at the time step n and n-1
            RHS(count_rhs) = u_n(i,j) - 0.5*Crank + 0.5*dt*(S(x(i),y(j),t) + S(x(i),y(j),t-dt)) + 0.5*(ghost_n + ghost_nm1);
        end
    end
end

% Finally, we are ready to build our sparse diffusion matrix A and solve
% for our new u_n
A = sparse(sparse_i,sparse_j,sparse_A);
sol = A\RHS;

for k = 1:RHS_numel
    u_n(untag_i(k),untag_j(k)) = sol(k);
end

end