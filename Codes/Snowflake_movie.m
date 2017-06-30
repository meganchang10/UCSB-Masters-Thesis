function [times_through_loop,phi,Frames] = Snowflake_movie(m,iterations,D,d0,epsilon,phi_eqn,tfinal)
tic
n = m;
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

contourf(x,y,phi',[-0.001 0.001],'k');
axis square
axis equal
colormap([0 0 0])
pause(dt);

dt = dx;
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

Frames(1) = getframe(gcf);
frames = 2;

times_through_loop = 0;
% tic
while t < tfinal

    count = count + 1;
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
%         toc
        count = 0;
        phi = Reinitialize(phi_np1,m,n,dx,dy,5*iterations_for_reinitialization);
%         tic
        dt
        times_through_loop
        t
    else
        phi = Reinitialize(phi_np1,m,n,dx,dy,iterations_for_reinitialization);
    end
    
    contourf(x,y,phi_np1',[-0.001 0.001],'k');
    axis square;
    axis equal;
    drawnow;

    Frames(frames) = getframe(gcf);
    pause(dt)
    frames = frames + 1;
    
    t = t+dt;
    
end

% video = VideoWriter('Snowflakeyolo .avi','Uncompressed AVI');
% open(video)
% writeVideo(video,Frames)
% close(video)

toc