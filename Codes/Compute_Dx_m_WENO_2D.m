function [Dx_m] = Compute_Dx_m_WENO_2D(u,dx,m,n,BCx)

% WENO is the Weighted Essentially Non-Oscillatory scheme

% The WENO schemes compute higher order one-sided first derivatives. For
% Dx_m, this is the side that would include i and i-1. WENO schemes choose
% between three possible stencils that can range from point i-3 to i+2 for
% Dx_m, selecting the stencil that will create the smoothest solution and
% avoid oscillatory results. This can involve only one stencil or a convex
% combination of all three stencils. To determine which case is best, we
% calculate the smoothness coefficients S1, S2, and S3, using the simple
% finite difference d1, d2, d3, d4, and d5. We then define the weight
% coefficients omega1, omega2, and omega3 which will determine how much
% each stencil u_x1, u_x2, and u_x3 will contribute to the final Dx_m
% solution.

Dx_m = zeros(m,n);
% Pad u with valid values that reflect the boundary conditions at the walls

u_padd = zeros(m+5,n);
for j = 1:n
    for i = 1:m
        u_padd(i+3,j) = u(i,j);
    end
    if (BCx==0)          % Padd u with Linear Extrapolation BC
        u_padd(  3,j) = 2*u_padd(  4,j) - u_padd(  5,j);
        u_padd(  2,j) = 2*u_padd(  3,j) - u_padd(  4,j);
        u_padd(  1,j) = 2*u_padd(  2,j) - u_padd(  3,j);
        u_padd(m+4,j) = 2*u_padd(m+3,j) - u_padd(m+2,j);
        u_padd(m+5,j) = 2*u_padd(m+4,j) - u_padd(m+3,j);
    elseif (BCx==1)      % Padd u with Periodic BC
        u_padd(  3,j) = u(m-1,j);
        u_padd(  2,j) = u(m-2,j);
        u_padd(  1,j) = u(m-3,j);
        u_padd(m+4,j) = u(2,j);
        u_padd(m+5,j) = u(3,j);
    end
end

% WENO formula
for j = 1:n
    for i = 4:m+3
        
        % Simple Finite Differences
        d1 = (u_padd(i-2,j) - u_padd(i-3,j))/dx;
        d2 = (u_padd(i-1,j) - u_padd(i-2,j))/dx;
        d3 = (u_padd(i  ,j) - u_padd(i-1,j))/dx;
        d4 = (u_padd(i+1,j) - u_padd(i  ,j))/dx;
        d5 = (u_padd(i+2,j) - u_padd(i+1,j))/dx;
        
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
        sum = alpha1 + alpha2 + alpha3;
        
        % Stencil Weights
        omega1 = alpha1/sum;
        omega2 = alpha2/sum;
        omega3 = 1-omega1-omega2;
        
        % Stencils
        u_x1 =  d1/3. -7*d2/6. +11*d3/6.;
        u_x2 = -d2/6. +5*d3/6. +d4/3.;
        u_x3 =  d3/3. +5*d4/6. -d5/6.;
        
        % Dx_m construction
        Dx_m(i-3,j) = omega1*u_x1 + omega2*u_x2 + omega3*u_x3;
    end
end
end