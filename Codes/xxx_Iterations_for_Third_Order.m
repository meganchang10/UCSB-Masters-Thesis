close all
clear all
clc

m = [81,161,321,641,1281];
%Plots will only use first 3 m's

exact_sol   = @(x,y) 0.5 - cos(x)*sin(y);


R = 0.35;
P = 0.28;
%Choose your level set shape here
% phi_initial = @(x,y,t) sqrt((x)^2 + (y)^2) - 0.5;                                                       %Circle
% phi_initial = @(x,y,t) abs(x) + abs(y) - 0.5;                                                           %Diamond
phi_initial = @(x,y,t) (x-0.02*sqrt(5))^2 + (y-0.02*sqrt(5))^2 - (0.5 + 0.2*sin(5*theta_maker(x,y)))^2; %Flower
% phi_initial = @(x,y,t) max(abs(x),abs(y)) - 0.4*sqrt(2);                                                %Square
% phi_initial = @(x,y,t) sqrt(x^2 + y^2) - (R + P*cos(9*theta_maker(x,y)));                             %Crazy Flower

%BC = 0 linear extrapolation
%BC = 1 periodic
BCx   = 1;
BCy   = 1;

xmin = -1; xmax = 1;
ymin = -1; ymax = 1;


iterations = [10 15 20 25 30 35 40 45 50 55 60 65];

error_inf = zeros(length(m),length(iterations));
error_1   = zeros(length(m),length(iterations));

for ii = 1:length(m)
for jj = 1:length(iterations)

tic
    n = m(ii);
    iteration = iterations(jj);
    disp(n)
    x = linspace(xmin,xmax,m(ii));
    y = linspace(ymin,ymax,n);
    dx = x(2)-x(1);
    dy = y(2)-y(1);
    
    dt = 0.5*min(dx,dy);
    
    phi   = zeros(n,n);
    T_out = zeros(n,n);
    exact = zeros(n,n);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:m(ii)
        for j = 1:n
            phi(i,j)   = phi_initial(x(i),y(j));
            exact(i,j) = exact_sol(x(i),y(j));
            if phi(i,j) <= 0
                T_out(i,j) = exact_sol(x(i),y(j));
            end
        end
    end
%     figure()
%     mesh(x,y,T_out)
%     title('Before')
    
    [nx,ny] = N_vector_generator(phi,n,n,dx,dy);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T_out_extrap =  Third_Order_Extrapolation(phi,T_out,n,n,dx,dy,iteration);
%     
%     figure()
%     mesh(x,y,T_out_extrap)
%     title('Extrapolated')
%     hold on
%     contour(x,y,phi,[-0.001,0.001])
    
    exact_plot = zeros(n,n);
    
    count = 0;
    for i = 1:m(ii)
        for j = 1:n
            exact_plot(i,j) = exact_sol(x(i),y(j));
            if abs(phi(i,j)) < 2*dx
                if phi(i,j) > 0
                    count = count + 1;
                end
            end
        end
    end
    
    error = zeros(count,1);
    count = 0;
    error_plot = zeros(n,n);
    for i = 1:m(ii)
        for j = 1:n
            error_plot(i,j) = exact_plot(i,j) - T_out_extrap(i,j);
            if abs(phi(i,j)) < 2*dx
                if phi(i,j) > 0
                    count = count + 1;
                    error(count) = exact(i,j) - T_out_extrap(i,j);
                   
                end
            end
        end
    end
    
    %THIS IS NORM INF ONLY FOR VECTORS
    error_1(ii,jj) = norm(error,1)*dx*dy;
    error_inf(ii,jj) = norm(error,Inf);
    
    %THIS IS NORM INF ONLY FOR MATRICES
    %error_inf(ii) = max(max(abs(error)));
    %error_1(ii) = sum(sum(abs(error), 1), 2)*dx*dy;
    
toc
end

end

x = iterations;
figure(1)
plot(x,log(error_inf(1,:)),'r-o',x,log(error_inf(2,:)),'b-o',x,log(error_inf(3,:)),'g-o',x,log(error_inf(4,:)),'k-o',x,log(error_inf(5,:)),'m-o')
legend('n = 81','n = 161','n = 321','n = 641','n = 1281')
xlabel('iterations')
ylabel('log(error)')
title('Error Analysis L^\infty')

figure(2)
plot(x,log(error_1(1,:)),'r-o',x,log(error_1(2,:)),'b-o',x,log(error_1(3,:)),'g-o',x,log(error_1(4,:)),'k-o',x,log(error_1(5,:)),'m-o')
legend('n = 81','n = 161','n = 321','n = 641','n = 1281')
xlabel('iterations')
ylabel('log(error)')
title('Error Analysis L^1')