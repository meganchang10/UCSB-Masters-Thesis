close all

m = [81,161,321,642];

exact_sol   = @(x,y) 0.5 - cos(x)*sin(y);

R = 0.35;
P = 0.28;
%Choose your level set shape here
phi_initial = @(x,y,t) sqrt((x)^2 + (y)^2) - 0.5;                                                       %Circle
% phi_initial = @(x,y,t) abs(x) + abs(y) - 0.5;                                                           %Diamond
% phi_initial = @(x,y,t) (x-0.02*sqrt(5))^2 + (y-0.02*sqrt(5))^2 - (0.5 + 0.2*sin(5*theta_maker(x,y)))^2; %Flower
% phi_initial = @(x,y,t) max(abs(x),abs(y)) - 0.4*sqrt(2);                                                %Square
% phi_initial = @(x,y,t) sqrt(x^2 + y^2) - (R + P*cos(9*theta_maker(x,y)));                             %Crazy Flower

R = 0.15;
P = 0.1;
% phi_initial = @(x,y,t) sqrt(x^2 + y^2) - (R + P*cos(5*theta_maker(x,y)));

%BC = 0 linear extrapolation
%BC = 1 periodic
BCx   = 0;
BCy   = 0;

xmin = -1; xmax = 1;
ymin = -1; ymax = 1;
error_inf = zeros(length(m),1);
error_1   = zeros(length(m),1);

n = m(1);
    x = linspace(xmin,xmax,n);
    y = linspace(ymin,ymax,n);
    dx = x(2)-x(1);
    dy = y(2)-y(1);
    
dist= 2*dx;

for ii = 1:length(m)
tic
    n = m(ii);
    disp(n)
    x = linspace(xmin,xmax,n);
    y = linspace(ymin,ymax,n);
    dx = x(2)-x(1);
    dy = y(2)-y(1);
    
    phi   = zeros(n,n);
    T_out = zeros(n,n);
    exact = zeros(n,n);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:n
        for i = 1:n
            phi(i,j)   = phi_initial(x(i),y(j));
            exact(i,j) = exact_sol(x(i),y(j));
            if phi(i,j) <= eps
                T_out(i,j) = exact_sol(x(i),y(j));
            end
        end
    end
    phi0 = phi;
    
    dt = 0.5*min(dx,dy);

%     phi = Reinitialize_new2_test(phi,n,n,dx,dy,dt,BCx,BCy,n/2);    
%     phi = Reinitialize_new_no(phi,n,n,dx,dy,dt,BCx,BCy,n/2);
%     phi = Reinitialize_new_no_WENO(phi,n,n,dx,dy,dt,n/2);
%     [phi] = Reinitialize_new(phi,n,n,dx,dy,dt,BCx,BCy,n/2);

    phi = Reinitialize_4_a_new_hope(phi,n,n,dx,dy,n/2);

    iterations = 50;
%     T_out_extrap = Third_Order_Extrapolation(phi,T_out,nx,ny,n,n,dx,dy,BCx,BCy,n);
    T_out_extrap = Third_Order_Extrapolation(phi,T_out,n,n,dx,dy,iterations);

    %Preallocate size of our error vector
    count = 0;
    for i = 1:n
        for j = 1:n% figure()
            if abs(phi(i,j)) < 2*dx
%              if abs(phi(i,j)) < dist
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
            if abs(phi(i,j)) < 2*dx
%              if abs(phi(i,j)) < dist
                if phi(i,j) > 0
                    count = count + 1;
                    error(count) = exact(i,j) - T_out_extrap(i,j);
                end
            end
        end
    end
    
    %THIS IS NORM FOR VECTORS
    error_inf(ii) = norm(error,Inf);
    error_1(ii) = norm(error,1)*dx*dy;

    %THIS IS NORM INF ONLY FOR MATRICES
    error_inf(ii) = max(max(abs(error)));
    error_1(ii) = sum(sum(abs(error), 1), 2)*dx*dy;
    
    error_plot = exact - T_out_extrap;
    
    figure()
    mesh(1:n,1:n,error_plot')
    hold on 
    contour(1:n,1:n,phi0',[-0.001,0.001],'r')
    contour(1:n,1:n,phi',[-0.001,0.001],'c')
    
    title('error plot')
    xlabel('x')
    ylabel('y')

toc
end

order(error_1,error_inf)

figure()
plot(log(m),log(error_inf),'r*')
xlabel('log(N)')
ylabel('log(error)')
title('Error Analysis L^\infty')

figure()
plot(log(m),log(error_1),'r*')
xlabel('log(N)')
ylabel('log(error)')
title('Error Analysis L^1')
