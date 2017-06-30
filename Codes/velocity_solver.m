function [u,v,Vn] = velocity_solver(phi,T_in,T_out,D_in,D_out,m,n,dx,dy)

% This function solves for the 'interface' velocity throughout the domain.
% However, we are only interested in the velocity that is actually at the
% interface.

% This velocity solver calculates second-order accurate first derivatives
% of the temperature gradient in both T_in and T_out which have been
% extrapolated in a band around the interface. We can then multiply these
% values by their respective diffusion coefficients D_in and D_out to find
% the magnitude of the velocity, and then multiply these by the normal
% vectors to determinethe velocity componenets u and v

% Below that, we have commented out another method to calculate the
% velocity components. Both are valid in their own way. The first one which
% we will use implements a physical condition that imposes the laws of
% physics. The second ensures greater numerical stability which is of use
% at higher grid resolutions.

% Vn is the magnitude of the interface velocity and u and v are the
% velocity components in the x and y direction, respectively. We have both
% as outputs because the interface velocity can be used in boundary
% condition calculations at the interface in order to create anisotropic
% effects.

u       = zeros(m,n);
v       = zeros(m,n);
Vn      = zeros(m,n);

[nx,ny] = N_vector_generator(phi,m,n,dx,dy);

for i = 2:m-1
    for j = 2:n-1
        DTinDx  = ( T_in(i+1,j  )  -T_in(i-1,j  ))/(2*dx);
        DTinDy  = ( T_in(i  ,j+1)  -T_in(i  ,j-1))/(2*dy);
        DToutDx = (T_out(i+1,j  ) -T_out(i-1,j  ))/(2*dx);
        DToutDy = (T_out(i  ,j+1) -T_out(i  ,j-1))/(2*dy);
        Vn(i,j) = -(D_out*DToutDx - D_in*DTinDx)*nx(i,j) - (D_out*DToutDy - D_in*DTinDy)*ny(i,j);
        u(i,j)  = Vn(i,j)*nx(i,j);
        v(i,j)  = Vn(i,j)*ny(i,j);
%         u(i,j) = - (D_out*DToutDx(i,j) - D_in*DTinDx(i,j));
%         v(i,j) = - (D_out*DToutDy(i,j) - D_in*DTinDy(i,j));
    end
end

end