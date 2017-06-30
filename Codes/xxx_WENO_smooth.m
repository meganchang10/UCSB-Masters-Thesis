%Testing_Weno_Accuracy
close all
clear all
clc


m = [81,161,321];
% m = [81,161,321];
% m = [41,81,161];
% m = 50;
%BC = 0 linear extrapolation
%BC = 1 periodic
BCx = 1;
BCy = 1;

phi_function = @(x,y) cos(x)*sin(y);
phi_dx = @(x,y) -sin(x)*sin(y);
phi_dy = @(x,y) cos(x)*cos(y);

error_inf_dxm = zeros(length(m),1);
error_1_dxm = zeros(length(m),1);
error_inf_dxp = zeros(length(m),1);
error_1_dxp = zeros(length(m),1);
error_inf_dym = zeros(length(m),1);
error_1_dym = zeros(length(m),1);
error_inf_dyp = zeros(length(m),1);
error_1_dyp = zeros(length(m),1);

for ii = 1:length(m)
    n = m(ii);
    phi = zeros(n,n);
    exact_dx = zeros(n,n);
    exact_dy = zeros(n,n);
    
    x = linspace(0,2*pi,n);
    y = linspace(0,2*pi,n);
    dx = x(2) - x(1);
    dy = y(2) - y(1);
    for i = 1:n
        for j = 1:n
            phi(i,j) = phi_function(x(i),y(j));
            exact_dx(i,j) = phi_dx(x(i),y(j));
            exact_dy(i,j) = phi_dy(x(i),y(j));
        end
    end

%     figure()
%     contour(x,y,phix)
%     title('x')
%     xlabel('x')
%     ylabel('y')
% 
%     figure()
%     mesh(x,y,phix)
%     title('x')
%     xlabel('x')
%     ylabel('y')
% 
%     figure()
%     contour(x,y,phiy)
%     title('not x')
%     xlabel('x')
%     ylabel('y')
% 
%     figure()
%     mesh(x,y,phiy)
%     title('not x')
%     xlabel('x')
%     ylabel('y')
    
    [Dx_m] = Compute_Dx_m_WENO_2D(phi,dx,n,n,BCx);
    [Dx_p] = Compute_Dx_p_WENO_2D(phi,dx,n,n,BCx);
    [Dy_m] = Compute_Dy_m_WENO_2D(phi,dy,n,n,BCy);
    [Dy_p] = Compute_Dy_p_WENO_2D(phi,dy,n,n,BCy);
    
    error_dxm = Dx_m-exact_dx;
    error_dxp = Dx_p-exact_dx;
    error_dym = Dy_m-exact_dy;
    error_dyp = Dy_p-exact_dy;
    
    error_inf_dxm(ii) = max(max(abs(error_dxm)));
    error_1_dxm(ii) = sum(sum(abs(error_dxm), 1), 2)*dx*dy;
    
    error_inf_dxp(ii) = max(max(abs(error_dxp)));
    error_1_dxp(ii) = sum(sum(abs(error_dxp), 1), 2)*dx*dy;
    
    error_inf_dym(ii) = max(max(abs(error_dym)));
    error_1_dym(ii) = sum(sum(abs(error_dym), 1), 2)*dx*dy;
    
    error_inf_dyp(ii) = max(max(abs(error_dyp)));
    error_1_dyp(ii) = sum(sum(abs(error_dyp), 1), 2)*dx*dy;
    
end

% figure()
% mesh(x,y,Dx_m_clean)
% title('Dxm')
% xlabel('x')
% ylabel('y')
% 
% figure()
% contour(x,y,Dx_p_clean)
% title('Dxp')
% xlabel('x')
% ylabel('y')
% 
% figure()
% mesh(x,y,Dy_m_clean)
% title('Dnotxm')
% xlabel('x')
% ylabel('y')
% 
% figure()
% contour(x,y,Dy_p_clean)
% title('Dnotxp')
% xlabel('x')
% ylabel('y')
% 

% 
% figure()
% plot(log(m),log(error_inf_dxm),'r*')
% xlabel('log(N)')
% ylabel('log(error)')
% title('Error Analsis L^\infty for Dxm')
% 
% figure()
% plot(log(m),log(error_1_dxm),'r*')
% xlabel('log(N)')
% ylabel('log(error)')
% title('Error Analsis L^1 for Dxm')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure()
% plot(log(m),log(error_inf_dxp),'r*')
% xlabel('log(N)')
% ylabel('log(error)')
% title('Error Analsis L^\infty for Dx+')
% 
% figure()
% plot(log(m),log(error_1_dxp),'r*')
% xlabel('log(N)')
% ylabel('log(error)')
% title('Error Analsis L^1 for Dxp')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure()
% plot(log(m),log(error_inf_dym),'r*')
% xlabel('log(N)')
% ylabel('log(error)')
% title('Error Analsis L^\infty for Dwhm')
% 
% figure()
% plot(log(m),log(error_1_dym),'r*')
% xlabel('log(N)')
% ylabel('log(error)')
% title('Error Analsis L^1 for Dwhm')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure()
% plot(log(m),log(error_inf_dyp),'r*')
% xlabel('log(N)')
% ylabel('log(error)')
% title('Error Analsis L^\infty for Dwh+')
% 
% figure()
% plot(log(m),log(error_1_dyp),'r*')
% xlabel('log(N)')
% ylabel('log(error)')
% title('Error Analsis L^1 for Dyp')

disp('Dxm Order')
order(error_1_dxm,error_inf_dxm)

disp('Dxp Order')
order(error_1_dxp,error_inf_dxp)

disp('Dym Order')
order(error_1_dym,error_inf_dym)

disp('Dyp Order')
order(error_1_dyp,error_inf_dyp)