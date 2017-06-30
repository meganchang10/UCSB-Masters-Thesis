close all

tic
m = [81,161];
m = [51,101];
m = 101;

%For testing just need to input exact_sol
%and calculate corresponding Source Term S
phi_initial   = @(x,y,t) (x-0.5)^2 + (y-0.75)^2 - 0.15^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tfinal = 1;
x0 = 0;     xf = 1;
y0 = 0;     yf = 1;
error_inf = zeros(length(m),1);
error_1   = zeros(length(m),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_steps = zeros(length(m),1);
for ii = 1:length(m)
    
    n = m(ii);
    disp(n)
    x = linspace(x0,xf,n);
    y = linspace(y0,yf,n);
    t = 0;
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    dt = dx;
    
    phi = zeros(n,n);
    u   = zeros(n,n);
    v   = zeros(n,n);
    for j = 1:n
        for i = 1:n
            phi(i,j) = phi_initial(x(i),y(j),t);
            u(i,j) = -sin(pi*x(i))^2*sin(2*pi*y(j));
            v(i,j) =  sin(pi*y(j))^2*sin(2*pi*x(i));
        end
    end
    
    phi_0 = phi;
    
    figure(1)
    contour(x,y,phi',[-0.0001,0.0001],'b')
    
    while t < tfinal
        if t + dt > tfinal
            dt = tfinal - t;
        end
        t = t + dt;
        phi = oneStepAdvection(phi,u,v,n,n,dx,dy,dt,0,0);
    end
    
    figure()
    contour(x,y,phi',[-0.0001,0.0001],'b')
    
    dt = dx;
    while t < 2*tfinal
        if t + dt > 2*tfinal
            dt = 2*tfinal - t;
        end
        t = t + dt;
        phi = oneStepAdvection(phi,-u,-v,n,n,dx,dy,dt,0,0);
    end
    
    figure(1)
    hold on
    contour(x,y,phi',[-0.0001,0.0001],'r')
    
    %THIS IS NORM INF ONLY FOR VECTORS
    %error_1(i) = norm(exact-u_n,1)*dx*dy;
    %error_inf(i) = norm(exact-u_n,Inf);
    %THIS IS NORM INF ONLY FOR MATRICES
    error_1(ii) = sum(sum(abs(phi_0-phi), 1), 2)*dx*dy;
    error_inf(ii) = max(max(abs(phi_0-phi)));

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