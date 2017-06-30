function [times_through_loop,phi] = Snowflake_fun(m,iterations,D,d0,epsilon,dt_plot,phi_eqn,tfinal,figure_number,time_steps_or_nah)

n = m;
%not sure how man iteraations I really need
iterations_for_extrapolation    = iterations;
iterations_for_reinitialization = iterations;
iterations_for_constant_extrap  = iterations;

D_in  =     D;
D_out =     D;
T_in   =    0;
T_out  = -0.5;

%BC = 0 linear extrapolation
%BC = 1 periodic
BCx   = 0;
BCy   = 0;
S = @(x,y,t) 0;

xmin = -1.5; xmax = 1.5;
ymin = -1.5; ymax = 1.5;
x = linspace(xmin,xmax,m);
y = linspace(ymin,ymax,n);
dx = x(2)-x(1);
dy = y(2)-y(1);

t = 0;
dt = 0;

phi     = zeros(m,n);
T_in_n  = zeros(m,n);
T_out_n = zeros(m,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:n
    for i = 1:m  
        phi(i,j) = phi_eqn(x(i),y(j));
        if phi(i,j) < 0
            T_in_n(i,j)  = T_in;
        end
        if phi(i,j) > 0
            T_out_n(i,j) = T_out;
        end
    end
end

phi = Reinitialize(phi,m,n,dx,dy,3*iterations_for_reinitialization);

figure(figure_number)
contour(x,y,phi',[-0.001 0.001]);
pause(dt);

dt = dx;
%BC_interface_np1 = T_in*ones(m,n);
%BC_interface = BC_interface_np1;

BC_interface_np1 = make_BC_interface(phi,d0,epsilon,m,n,dx,dy);
BC_interface = BC_interface_np1;
for Queen_Elsa = 1:iterations*2
    T_in_n  = oneStepDiffusion(T_in_n , phi,m,n,D_in ,S,BC_interface_np1,T_in,BC_interface ,x,y,dx,dy,dt,t);
    T_out_n = oneStepDiffusion(T_out_n,-phi,m,n,D_out,S,BC_interface_np1,T_out,BC_interface,x,y,dx,dy,dt,t);
end

T_in_np1 = T_in_n;
T_out_np1 = T_out_n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 0;
count_plot = 0;
times_through_loop = 0;
t_plot = dt_plot;

tic
while t < tfinal

    count = count + 1;
    count_plot = count_plot + 1;
    times_through_loop = times_through_loop + 1;

    T_in_n  = Third_Order_Extrapolation( phi,T_in_np1 ,m,n,dx,dy,iterations_for_extrapolation);
    T_out_n = Third_Order_Extrapolation(-phi,T_out_np1,m,n,dx,dy,iterations_for_extrapolation);
    
    [u,v] = velocity_solver(phi,T_in_n,T_out_n,D_in,D_out,m,n,dx,dy);
    
    u = Constant_Extrapolation(u, phi,m,n,dx,dy,iterations_for_constant_extrap);
    u = Constant_Extrapolation(u,-phi,m,n,dx,dy,iterations_for_constant_extrap);
    v = Constant_Extrapolation(v, phi,m,n,dx,dy,iterations_for_constant_extrap);
    v = Constant_Extrapolation(v,-phi,m,n,dx,dy,iterations_for_constant_extrap);

    Umax = max(max(sqrt(u.^2 + v.^2)));
    dt = 0.5*min(dx,dy)/Umax;
    phi_np1 = oneStepAdvection(phi,u,v,m,n,dx,dy,dt,BCx,BCy);

    BC_interface = BC_interface_np1;
    BC_interface_np1 = make_BC_interface(phi_np1,d0,epsilon,m,n,dx,dy);

    T_in_np1  = oneStepDiffusion(T_in_n , phi_np1,m,n,D_in ,S,BC_interface_np1,T_in,BC_interface ,x,y,dx,dy,dt,t);
    T_out_np1 = oneStepDiffusion(T_out_n,-phi_np1,m,n,D_out,S,BC_interface_np1,T_out,BC_interface,x,y,dx,dy,dt,t);
   
    if count == 5
        toc
        count = 0;
        phi = Reinitialize(phi_np1,m,n,dx,dy,5*iterations_for_reinitialization);
        tic
        dt
        times_through_loop
    else
        phi = Reinitialize(phi_np1,m,n,dx,dy,iterations_for_reinitialization);
    end
    
    if t > t_plot 
        t_plot = t_plot + dt_plot;
        figure(figure_number)
        hold on
        contour(x,y,phi_np1',[-0.001 0.001]);
        title(['t = ' num2str(t+dt)])
        axis square;
        axis equal;
        pause(dt);
        if time_steps_or_nah == 1            
            figure()
            contour(x,y,phi_np1',[-0.001 0.001],'k');
            title(['t = ' num2str(t+dt)])
            axis square;
            axis equal;
            pause(dt);
        end
    end

    t = t+dt;
    
end
toc