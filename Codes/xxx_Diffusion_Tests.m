close all

tic
m = [81,161,321];

D = 2; %Diffusion Coefficient

%For testing just need to input exact_sol
%and calculate corresponding Source Term S
exact_sol   = @(x,y,t) cos(x)*sin(x)*exp(y)*sin(t);
S           = @(x,y,t) cos(x)*sin(x)*exp(y)*cos(t)-D*(-2*sin(2*x)*exp(y)*sin(t)+cos(x)*sin(x)*exp(y)*sin(t)); %Source Function

R = 0.35;
P = 0.28;
%Choose your level set shape here
phi_initial = @(x,y,t) sqrt((x)^2 + (y)^2) - 0.5;                                                       %Circle
% phi_initial = @(x,y,t) abs(x) + abs(y) - 0.5;                                                           %Diamond
% phi_initial = @(x,y,t) (x-0.02*sqrt(5))^2 + (y-0.02*sqrt(5))^2 - (0.5 + 0.2*sin(5*theta_maker(x,y)))^2; %Flower
% phi_initial = @(x,y,t) max(abs(x),abs(y)) - 0.4*sqrt(2);                                                %Square
% phi_initial = @(x,y,t) sqrt(x^2 + y^2) - (R + P*cos(9*theta_maker(x,y)));                             %Crazy Flower

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BC = exact_sol;
Initial_Data = exact_sol;

tfinal = 0.1;
x0 = -1;     xf = 1;
y0 = -1;     yf = 1;
error_inf = zeros(length(m),1);
error_1   = zeros(length(m),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:length(m)
    
    n = m(ii);
    x = linspace(x0,xf,n);
    y = linspace(y0,yf,n);
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    t = 0;
    
    BC_interface = zeros(n,n);
    BC_interface_old = zeros(n,n);
    BC_domain = zeros(n,n);
    
    dt = dx;
    
    phi = zeros(n,n);
    u_n = zeros(n,n);
    for i = 1:n
        for j = 1:n
            phi(i,j) = phi_initial(x(i),y(j),0);
            if phi(i,j) < 0
                u_n(i,j) = Initial_Data(x(i),y(j),0);
            end
        end
    end
    
    count = 0;
    while t < tfinal
        if t + dt > tfinal
            dt = tfinal - t;
        end
        count = count + 1;
        t = t + dt;
        for j = 1:n
            for i = 1:n
                BC_interface(i,j) = exact_sol(x(i),y(j),t);
                BC_interface_old(i,j) = exact_sol(x(i),y(j),t-dt);
            end
        end
        BC_domain = BC_interface;
        u_n = oneStepDiffusion_BC_old(u_n,phi,n,n,D,S,BC_interface,BC_domain,BC_interface_old,x,y,dx,dy,dt,t);
%         u_n = oneStepDiffusion_analytical(u_n,phi,n,n,D,S,BC,BC,x,y,dx,dy,dt,t);
    end
    
    exact = analyticalsolution_diffusion(n,x0,xf,y0,yf,tfinal,exact_sol,phi_initial);
    
    %THIS IS NORM INF ONLY FOR VECTORS
    %error_inf(i) = norm(exact-u_n,Inf);
    %error_1(i) = norm(exact-u_n,1)*dx*dy;

    %THIS IS NORM INF ONLY FOR MATRICES
    error_inf(ii) = max(max(abs(exact-u_n)));
    error_1(ii) = sum(sum(abs(exact-u_n), 1), 2)*dx*dy;
    
    %3D Plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure()
    xx = zeros(n,n);
    yy = zeros(n,n);
    u_n_high = zeros(n,n);
    for i = 1:n
        for j = 1:n
            xx(:,j) = x(j);
            yy(i,:) = y(i);
            if not(u_n(i,j) == 0)
                u_n_high(i,j) = u_n(i,j) + 0.2;
            end
        end
    end
    plot3(xx,yy,u_n_high,'ro')
    xlabel('x')
    ylabel('y')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Mesh Plot
    figure()
    mesh(x,y,u_n)
    xlabel('x')
    ylabel('y')
    
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
toc