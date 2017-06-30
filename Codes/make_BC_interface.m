function [BC_interface] = make_BC_interface(phi,d0,epsilon,m,n,dx,dy)
% epsilon = anisotropy strength
% alpha   = angle b/n normal at interface and x axis
% d0      = some coefficient
% k       = curvature

BC_interface = zeros(m,n);
[nx,ny] = N_vector_generator(phi,m,n,dx,dy);

for j = 2:n-1
    for i = 2:m-1
        alpha     = atan(ny(i,j)/nx(i,j));
        epsilon_c = d0*(1-15*epsilon*cos(4*alpha));
        
        n_x  = (nx(i+1,j  ) - nx(i-1,j  ))/(2*dx);
        n_y  = (ny(i  ,j+1) - ny(i  ,j-1))/(2*dy);
        curvature = n_x + n_y;
        
        BC_interface(i,j) = -epsilon_c*curvature;
        if BC_interface(i,j) < 0
            BC_interface(i,j) = 0;
        end
    end
end

end