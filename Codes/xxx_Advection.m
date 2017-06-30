close all
clear all
clc

m = [81,161];

% Boundary conditions:
% BC=0 for linear extrapolation
% BC=1 for periodic

BCx = 1;
BCy = 1;

tfinal = 1;
cfl    = 0.5;

error_inf = zeros(1,length(m));
error_1 = zeros(1,length(m));

for ii = 1:length(m)
    
    tic
    n = m(ii);
    disp(n);
    x = linspace(-pi,pi,n);
    y = linspace(-pi,pi,n);
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    dt = cfl*dx;
    
    phi        = zeros(n,n);
    exact_sol  = zeros(n,n);
    u          = zeros(n,n);
    v          = zeros(n,n);
    error      = zeros(n,n);
    
    for i = 1:n
        for j = 1:n
            u(i,j) = 1;
            v(i,j) = 1;
        end
    end
    
    xc = 0;
    yc = 0;
    xf = 1;
    yf = 1;
    Radius = 1;
    
    for j = 1:n
        for i = 1:n
            phi(i,j)       = cos(x(i)-xc)*sin(y(j)-yc);
            exact_sol(i,j) = cos(x(i)-xf)*sin(y(j)-yf);
        end
    end
    
    t = 0;
    
    % Evolve the solution in time
    while t < tfinal
        if (t + dt > tfinal)
            dt = tfinal - t;
        end
        t = t + dt;
        phi = oneStepAdvection(phi,u,v,n,n,dx,dy,dt,BCx,BCy);
    end


    for i = 1:n
        for j = 1:n
            error(i,j) = phi(i,j) - exact_sol(i,j);
        end
    end

    %THIS IS NORM INF ONLY FOR MATRICES
    error_1(ii) = sum(sum(abs(error), 1), 2)*dx*dy;
    error_inf(ii) = max(max(abs(error)));
    
    %THIS IS NORM 1 ONLY FOR VECTORS
    %error_1(ii) = norm(error,1)*dx*dy;
    %error_inf(ii) = norm(error,Inf);
    
    toc
end

order(error_1,error_inf)

figure()
plot(log(m),log(error_inf),'r*')
xlabel('log(N)')
ylabel('log(error)')
title('Error Analsis L^\infty')

figure()
plot(log(m),log(error_1),'r*')
xlabel('log(N)')
ylabel('log(error)')
title('Error Analsis L^1')
