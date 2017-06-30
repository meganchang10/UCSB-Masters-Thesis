function phi_np1 = oneStepAdvection(phi_n, u, v, m, n, dx, dy, dt, BCx, BCy)

% Uses WENO schemes to calculate higher-order one-sided finite differences.
% Then implements a TVD RK3 method which uses 3 Euler steps and 2 averaging
% steps to create a third-order accurate solution for our advection.
% Depending on the direction of the propagation of information (i.e. if 
% u > 0 or u < 0, we will use Dx_m or Dx_p, respectively. This function
% takes as an input the x and y velocity components u and v which are
% calculated at the interface. This interface velocity is then used to move
% the entire level set. However, it is really only the zero level set that
% we care about.

phi_np1     = zeros(m,n);
phi_np2     = zeros(m,n);
phi_np3half = zeros(m,n);

%% Implement a TVD RK3:
Dx_m = Compute_Dx_m_WENO_2D(phi_n,dx,m,n,BCx);
Dx_p = Compute_Dx_p_WENO_2D(phi_n,dx,m,n,BCx);
Dy_m = Compute_Dy_m_WENO_2D(phi_n,dy,m,n,BCy);
Dy_p = Compute_Dy_p_WENO_2D(phi_n,dy,m,n,BCy);

% First Euler Step to find phi_np1
for j = 1:n
    for i = 1:m
        if(u(i,j)>0)
            phi_np1(i,j) = phi_n(i,j)-dt*u(i,j)*Dx_m(i,j);
        else
            phi_np1(i,j) = phi_n(i,j)-dt*u(i,j)*Dx_p(i,j);
        end
    end
end

for j = 1:n
    for i = 1:m
        if(v(i,j)>0)
            phi_np1(i,j) = phi_np1(i,j)-dt*v(i,j)*Dy_m(i,j);
        else
            phi_np1(i,j) = phi_np1(i,j)-dt*v(i,j)*Dy_p(i,j);
        end
    end
end

Dx_m = Compute_Dx_m_WENO_2D(phi_np1,dx,m,n,BCx);
Dx_p = Compute_Dx_p_WENO_2D(phi_np1,dx,m,n,BCx);
Dy_m = Compute_Dy_m_WENO_2D(phi_np1,dy,m,n,BCy);
Dy_p = Compute_Dy_p_WENO_2D(phi_np1,dy,m,n,BCy);

% Second Euler Step to find phi_np2
for j = 1:n
    for i = 1:m
        if(u(i,j)>0)
            phi_np2(i,j) = phi_np1(i,j)-dt*u(i,j)*Dx_m(i,j);
        else
            phi_np2(i,j) = phi_np1(i,j)-dt*u(i,j)*Dx_p(i,j);
        end
    end
end

for j = 1:n
    for i = 1:m
        if(v(i,j)>0)
            phi_np2(i,j) = phi_np2(i,j)-dt*v(i,j)*Dy_m(i,j);
        else
            phi_np2(i,j) = phi_np2(i,j)-dt*v(i,j)*Dy_p(i,j);
        end
    end
end

% Average to get unp1half:
phi_np1half=.75*phi_n + .25*phi_np2;

Dx_m = Compute_Dx_m_WENO_2D(phi_np1half,dx,m,n,BCx);
Dx_p = Compute_Dx_p_WENO_2D(phi_np1half,dx,m,n,BCx);
Dy_m = Compute_Dy_m_WENO_2D(phi_np1half,dy,m,n,BCy);
Dy_p = Compute_Dy_p_WENO_2D(phi_np1half,dy,m,n,BCy);

% Third Euler Step to go from unp1half to unp3half
for j = 1:n
    for i = 1:m
        if(u(i,j)>0)
            phi_np3half(i,j) = phi_np1half(i,j)-dt*u(i,j)*Dx_m(i,j);
        else
            phi_np3half(i,j) = phi_np1half(i,j)-dt*u(i,j)*Dx_p(i,j);
        end
    end
end


for j = 1:n
    for i = 1:m
        if(v(i,j)>0)
            phi_np3half(i,j) = phi_np3half(i,j)-dt*v(i,j)*Dy_m(i,j);
        else
            phi_np3half(i,j) = phi_np3half(i,j)-dt*v(i,j)*Dy_p(i,j);
        end
    end
end

% Average to get unp1:
phi_np1=1./3*phi_n + 2./3*phi_np3half;
end