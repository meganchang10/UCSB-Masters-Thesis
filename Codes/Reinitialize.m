function phi = Reinitialize(phi,m,n,dx,dy,iterations)

dt = 0.5*min(dx,dy);
phi_np1 = zeros(m,n);
phi_np2 = zeros(m,n);


thetax_p = zeros(m,n);
thetax_m = zeros(m,n);
thetay_p = zeros(m,n);
thetay_m = zeros(m,n);
dtau     = dt*ones(m,n);

% Find thetas for entire grid and set them static for the rest of the
% iterations, no need to recalculate, anymore

% Reinitialization makes gradient linear to remove curvature, but must not
% move the zero level set which keeps track of the interface. We can
% calculate the position of the interface by computing theta, the
% fractional distance of the interface (phi = 0) from the grid node. We
% will save these values outside of the iterative loop so that they do not
% change, and in this way we will maintain the zero level set. 

% Notice that one grid point can have up to 4 different theta values if the 
% interface were surrounding it. For instance if the interface passed
% between a point phi(i,j) and phi(i+1,j), we would keep track of this
% theta in thetax_p in the position thetax_p(i,j). Theta must always be
% positive.

% Additionally, here, we will reset phi to 0 if theta is too small since
% this will cause spikes in our calculations since theta is in the
% denominator and fractions will then leap to infinity

% Finally, we must also keep track of the new time step that we will
% use at each grid point to ensure stability. This stability is maintained
% by using dtau = 0.5*min(dx,dy,theta)

for j = 2:n-1
    for i = 2:m-1
        if phi(i,j)*phi(i+1,j) < 0
            theta = phi(i,j)/(phi(i,j)-phi(i+1,j));
            %If theta is tooooo small then resets phi to 0
            if theta < dx
                phi(i,j) = 0;
            else
                thetax_p(i,j) = theta*dx;
                dtau(i,j) = min(dtau(i,j),0.5*thetax_p(i,j));
            end
        end
        if phi(i,j)*phi(i-1,j) < 0
            theta = phi(i,j)/(phi(i,j)-phi(i-1,j));
            if theta < dx
                phi(i,j) = 0;
            else
                thetax_m(i,j) = theta*dx;
                dtau(i,j) = min(dtau(i,j),0.5*thetax_m(i,j));
            end
        end
        if phi(i,j)*phi(i,j+1) < 0
            theta = phi(i,j)/(phi(i,j)-phi(i,j+1));
            if theta < dy
                phi(i,j) = 0;
            else
                thetay_p(i,j) = theta*dy;
                dtau(i,j) = min(dtau(i,j),0.5*thetay_p(i,j));
            end
        end
        if phi(i,j)*phi(i,j-1) < 0
            theta = phi(i,j)/(phi(i,j)-phi(i,j-1));
            if theta < dy
                phi(i,j) = 0;
            else
                thetay_m(i,j) = theta*dy;
                dtau(i,j) = min(dtau(i,j),0.5*thetay_m(i,j));
            end
        end
    end
end

% Now, we enter our iterative process where we calculate derivative in both
% upwind directions for both x and y where Dx_p is the upwind direction
% containing points at i+1 and i and Dx_m contains points i-1 and i.
% Special care is taken at the points where i = 1 since Dx_m cannot be
% calculated because the point 0 would be needed and this is not a real
% point. Thus, we simply set Dx_m = 0 for these points. The same is done
% for i = m since Dx_p cannot be calcualated since the point m+1 doesn't
% exist. The same is done in the y directions for j = 1 and j = n.

% Then we apply the Godunov Hamiltonian equation as shown below to
% calculate phi_np1 for our first Euler Step. The same process is done to
% calculate the phi_np2 for our second Euler Step. These 2 values are then
% averaged to get a more accurate value for phi_np1. This is Runge Kutta

for t = 1:iterations    
    for j = 1:n
        for i = 1:m

            if i==1
                Dx_p = (phi(i+1,j)-phi(i,j))/dx;
                Dx_m = 0;
            elseif i==m
                Dx_p = 0;
                Dx_m = (phi(i,j)-phi(i-1,j))/dx;
            else
                if phi(i,j)*phi(i+1,j) < 0
                    Dx_p = -phi(i,j)/thetax_p(i,j);
                else
                    Dx_p = (phi(i+1,j)-phi(i,j))/dx;
                end
                if phi(i,j)*phi(i-1,j) < 0
                    Dx_m = phi(i,j)/thetax_m(i,j);
                else
                    Dx_m = (phi(i,j)-phi(i-1,j))/dx;
                end
            end
            
            if j==1
                Dy_p = (phi(i,j+1)-phi(i,j))/dy;
                Dy_m = 0;
            elseif j==n
                Dy_p = 0;
                Dy_m = (phi(i,j)-phi(i,j-1))/dy;
            else
                if phi(i,j)*phi(i,j+1) < 0
                    Dy_p = -phi(i,j)/thetay_p(i,j);
                else
                    Dy_p = (phi(i,j+1)-phi(i,j))/dy;
                end
                if phi(i,j)*phi(i,j-1) < 0
                    Dy_m = phi(i,j)/thetay_m(i,j);
                else
                    Dy_m = (phi(i,j)-phi(i,j-1))/dy;
                end
            end
            
            if phi(i,j) < 0
                a_p = max(Dx_p,0);
                b_m = min(Dx_m,0);
                c_p = max(Dy_p,0);
                d_m = min(Dy_m,0);
                H_g = sqrt(max(a_p^2,b_m^2) + max(c_p^2,d_m^2));
                phi_np1(i,j) = phi(i,j) + dtau(i,j)*(H_g - 1);
            elseif phi(i,j) > 0
                a_m = min(Dx_p,0);
                b_p = max(Dx_m,0); 
                c_m = min(Dy_p,0); 
                d_p = max(Dy_m,0); 
                H_g = sqrt(max(a_m^2,b_p^2) + max(c_m^2,d_p^2));
                phi_np1(i,j) = phi(i,j) - dtau(i,j)*(H_g - 1);
            end
        end
    end

    for j = 1:n
        for i = 1:m
            if i==1
                Dx_p = (phi_np1(i+1,j)-phi_np1(i,j))/dx;
                Dx_m = 0;
            elseif i==m
                Dx_p = 0;
                Dx_m = (phi_np1(i,j)-phi_np1(i-1,j))/dx;
            else
                if phi_np1(i,j)*phi_np1(i+1,j) < 0
                    Dx_p = -phi_np1(i,j)/thetax_p(i,j);
                else
                    Dx_p = (phi_np1(i+1,j)-phi_np1(i,j))/dx;
                end
                if phi_np1(i,j)*phi_np1(i-1,j) < 0
                    Dx_m = phi_np1(i,j)/thetax_m(i,j);
                else
                    Dx_m = (phi_np1(i,j)-phi_np1(i-1,j))/dx;
                end
            end
            
            if j==1
                Dy_p = (phi_np1(i,j+1)-phi_np1(i,j))/dy;
                Dy_m = 0;
            elseif j==n
                Dy_p = 0;
                Dy_m = (phi_np1(i,j)-phi_np1(i,j-1))/dy;
            else
                if phi_np1(i,j)*phi_np1(i,j+1) < 0
                    Dy_p = -phi_np1(i,j)/thetay_p(i,j);
                else
                    Dy_p = (phi_np1(i,j+1)-phi_np1(i,j))/dy;
                end
                if phi_np1(i,j)*phi_np1(i,j-1) < 0
                    Dy_m = phi_np1(i,j)/thetay_m(i,j);
                else
                    Dy_m = (phi_np1(i,j)-phi_np1(i,j-1))/dy;
                end
            end

            if phi_np1(i,j) < 0
                a_p = max(Dx_p,0);
                b_m = min(Dx_m,0);
                c_p = max(Dy_p,0);
                d_m = min(Dy_m,0);
                H_g = sqrt(max(a_p^2,b_m^2) + max(c_p^2,d_m^2));
                phi_np2(i,j) = phi_np1(i,j) + dtau(i,j)*(H_g - 1);
            elseif phi_np1(i,j) > 0
                a_m = min(Dx_p,0);
                b_p = max(Dx_m,0);
                c_m = min(Dy_p,0);
                d_p = max(Dy_m,0);
                H_g = sqrt(max(a_m^2,b_p^2) + max(c_m^2,d_p^2));
                phi_np2(i,j) = phi_np1(i,j) - dtau(i,j)*(H_g - 1);
            end
        end
    end

    phi = (phi + phi_np2)/2;
    
end


end