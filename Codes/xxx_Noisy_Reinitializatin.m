close all
clear all

% m = [81,161,321];
m = 61;


% Boundary conditions:
% BC=0 for linear extrapolation
% BC=1 for periodic
BCx=0;
BCy=0;

R = 0.1;
P = 0.5;
%Choose your level set shape here
phi_initial = @(x,y) (x)^2 + (y)^2 - 0.65;                                                       %Circle
% phi_initial = @(x,y) abs(x) + abs(y) - 0.5;                                                           %Diamond
% phi_initial = @(x,y) (x-0.02*sqrt(5))^2 + (y-0.02*sqrt(5))^2 - (0.5 + 0.2*sin(5*theta_maker(x,y)))^2; %Flower
% phi_initial = @(x,y) max(abs(x),abs(y)) - 0.5*sqrt(2);                                                %Square
% phi_initial = @(x,y) sqrt(x^2 + y^2) - (R + P*cos(5*theta_maker(x,y)));


error_1   = zeros(length(m),1);
error_inf = zeros(length(m),1);

for ii = 1:length(m)
    n = m(ii);
    
    x = linspace(-1,1,n);
    y = linspace(-1,1,n);
    dx = x(2)-x(1);
    dy = y(2)-y(1);
    dt = 0.5*min(dx,dy);
    
    phi = zeros(n,n);
    for i = 1:n
        for j = 1:n
            phi(i,j) = phi_initial(x(i),y(j)) + 0.9*rand*dx;
        end
    end
    
    phi_0 = phi;
    
    figure()
    mesh(x,y,phi');
    hold on
    contour(x,y,phi',[-0.001,0.001],'r')
    axis([-1 1 -1 1 -0.8 1])
    
    iterations = 200;
    
    phi = Reinitialize_4_a_new_hope(phi,n,n,dx,dy,iterations);
    
    figure()
    mesh(x,y,phi');
    hold on
    contour(x,y,phi',[-0.0001,0.0001],'r')
    axis([-1 1 -1 1 -0.8 0.6])
%     
%     figure()
%     hold on
%     contour(1:n,1:n,phi_0',[-0.0000001,0.0000001],'r')
%     contour(1:n,1:n,phi',[-0.0000001,0.0000001],'b')
%     xlabel('x')
%     ylabel('y')
%     for k = 1:n
%         plot(1:n,k*ones(1,n),'-k');
%         plot(k*ones(1,n),1:n,'-k')
%     end
%     
end