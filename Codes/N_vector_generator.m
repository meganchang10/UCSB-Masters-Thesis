function [nx,ny] = N_vector_generator(phi,m,n,dx,dy)

% Our normal vector calculates the vector that is perpendicular to the
% surface of phi at any given point.

nx    = zeros(m,n);
ny    = zeros(m,n);
for j = 2:n-1
    for i = 2:m-1
        phi_x     = (phi(i+1,j  ) - phi(i-1,j  ))/(2*dx);
        phi_y     = (phi(i  ,j+1) - phi(i  ,j-1))/(2*dy);
        norm      = sqrt(phi_x^2 + phi_y^2);
        nx(i,j)   = phi_x/norm;
        ny(i,j)   = phi_y/norm;
    end
end
end