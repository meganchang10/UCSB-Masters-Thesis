function [Dy_m]=Compute_Dy_m_WENO_2D(u,dy,m,n,BCy)

% WENO is the Weighted Essentially Non-Oscillatory scheme

% The WENO schemes compute higher order one-sided first derivatives. For
% Dy_m, this is the side that would include j and j-1. WENO schemes choose
% between three possible stencils that can range from point j-3 to j+2 for
% Dy_m, selecting the stencil that will create the smoothest solution and
% avoid oscillatory results. This can involve only one stencil or a convex
% combination of all three stencils. To determine which case is best, we
% calculate the smoothness coefficients S1, S2, and S3, using the simple
% finite difference d1, d2, d3, d4, and d5. We then define the weight
% coefficients alpha1, alpha2, and alpha3 which will determine how much
% each stencil u_y1, u_y2, and u_y3 will contribute to the final Dy_m
% solution.

Dy_m=zeros(m,n);
% Pad u with valid values that reflect the boundary conditions at the walls

u_padd=zeros(m,n+5);
for i = 1:m
    for j = 1:n
        u_padd(i,j+3) = u(i,j);
    end
    if (BCy==0)          % Padd u with Linear Extrapolation BC
        u_padd(i,  3) = 2*u_padd(i,  4) - u_padd(i,  5);
        u_padd(i,  2) = 2*u_padd(i,  3) - u_padd(i,  4);
        u_padd(i,  1) = 2*u_padd(i,  2) - u_padd(i,  3);
        u_padd(i,n+4) = 2*u_padd(i,n+3) - u_padd(i,n+2);
        u_padd(i,n+5) = 2*u_padd(i,n+4) - u_padd(i,n+3);
    elseif (BCy==1)      % Padd u with Periodic BC
        u_padd(i,  3) = u(i,m-1);
        u_padd(i,  2) = u(i,m-2);
        u_padd(i,  1) = u(i,m-3);
        u_padd(i,m+4) = u(i,2  );
        u_padd(i,m+5) = u(i,3  );
    end
end

% WENO formula
for i = 1:m
    for j = 4:n+3
        
        % Simple Finite Differences
        d1 = (u_padd(i,j-2) - u_padd(i,j-3))/dy;
        d2 = (u_padd(i,j-1) - u_padd(i,j-2))/dy;
        d3 = (u_padd(i,j  ) - u_padd(i,j-1))/dy;
        d4 = (u_padd(i,j+1) - u_padd(i,j  ))/dy;
        d5 = (u_padd(i,j+2) - u_padd(i,j+1))/dy;
        
        % Smoothness Coefficients 
        S1 = 13./12*(d1-2*d2+d3)^2 + .25*(d1-4*d2+3*d3)^2;
        S2 = 13./12*(d2-2*d3+d4)^2 + .25*(d2-d4)^2;
        S3 = 13./12*(d3-2*d4+d5)^2 + .25*(3*d3-4*d4+d5)^2;
        
        Dsquare = [d1^2 d2^2 d3^2 d4^2 d5^2];
        epsilon = 1e-6*max(Dsquare) + 1e-99;
        
        alpha1 = .1/((S1+epsilon)^2);
        alpha2 = .6/((S2+epsilon)^2);
        alpha3 = .3/((S3+epsilon)^2);
        
        % Convex Combination
        sum    = alpha1 + alpha2 + alpha3;
        
        % Stencil Weights
        omega1 = alpha1/sum;
        omega2 = alpha2/sum;
        omega3 = 1-omega1-omega2;
        
        % Stencils
        u_y1 =  d1/3. - 7*d2/6. + 11*d3/6.;
        u_y2 = -d2/6. + 5*d3/6. + d4/3.;
        u_y3 =  d3/3. + 5*d4/6. - d5/6.;
        
        % Dy_m construction
        Dy_m(i,j-3) = omega1*u_y1 + omega2*u_y2 + omega3*u_y3;
    end
end
end