close all
clear all
clc

m = [51,101,201,401];
m = 61;

error_inf = zeros(length(m),1);
error_1   = zeros(length(m),1);

% exact_sol = @(x,y) exp(cos(y)+x*sinh(y))*cos(sin(x));
exact_sol = @(x,y) 1;

for ii = 1:length(m)
    n = m(ii);
    
    x = linspace(-1,1,n); dx = x(2)-x(1);
    y = linspace(-1,1,n); dy = y(2)-y(1);
    dt = .5*min(dx, dy);
    
    phi   = zeros(n,n);
    Q     = zeros(n,n);
    
    % Boundary conditions:
    % BC=0 for linear extrapolation
    % BC=1 for periodic
    BCx=1;
    BCy=1;
    
    R = 0.5;
    exact = zeros(n,n);
    for i = 1:n
        for j = 1:n
            phi(i,j) = sqrt(x(i)^2 + y(j)^2) - R;
            if phi(i,j) < 0
                Q(i,j) = exact_sol(x(i),y(j));
                exact(i,j) = exact_sol(x(i),y(j));
            else
                Q(i,j) = 0;
                theta = theta_maker(x(i),y(j));
                exact(i,j) = exact_sol(R*cos(theta),R*sin(theta));
            end
        end
    end

    
    figure()
    hold on
    contour(x,y,phi',[-.001 .001],'r')
    mesh(x,y,Q');
    
    [nx,ny] = N_vector_generator(phi,n,n,dx,dy);
    
    iterations = 15;
    tic
    Q1 = Constant_Extrapolation(Q,phi,n,n,dx,dy,iterations);
    toc
    
    figure()
    hold on
    contour(x,y,phi',[-.001 .001],'r')
    mesh(x,y,Q1');
    
    error_plot = zeros(n,n);
    for j = 5:n-4
        for i = 5:n-4
            error_plot(i,j) = exact(i,j) - Q1(i,j);
        end
    end
    
    count = 0;
    for i = 1:n
        for j = 1:n% figure()
            if abs(phi(i,j)) < 4*dx*2^ii
                if phi(i,j) > 0
                    count = count + 1;
                end
            end
        end
    end
    
    error = zeros(count,1);
    count = 0;
    for i = 1:n
        for j = 1:n
            if abs(phi(i,j)) < 4*dx*2^ii
                if phi(i,j) > 0
                    count = count + 1;
                    error(count) = exact(i,j) - Q1(i,j);
                end
            end
        end
    end
    
    %THIS IS NORM FOR VECTORS
    error_inf(ii) = norm(error,Inf);
    error_1(ii) = norm(error,1)*dx*dy;
    
%     %THIS IS NORM INF ONLY FOR MATRICES
%     error_inf(ii) = max(max(abs(error)));
%     error_1(ii) = sum(sum(abs(error), 1), 2)*dx*dy;

figure()
hold on
contour(x,y,phi',[-.001 .001],'r')
mesh(x,y,error_plot');
title('Error')
    
end

order(error_1,error_inf)
