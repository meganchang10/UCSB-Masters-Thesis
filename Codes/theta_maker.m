function theta = theta_maker(x,y)

% Takes as an input x and y and returns angle

if x > 0 && y >= 0        %Quadrant 1
    theta = atan(y/x);
elseif x < 0 && y >= 0    %Quadrant 2
    theta  = pi - atan(abs(y/x));
elseif x < 0 && y <= 0    %Quadrant 3
    theta  = pi + atan(abs(y/x));
elseif x > 0 && y <= 0    %Quadrant 4
    theta  = 2*pi - atan(abs(y/x));    
elseif x == 0 && y >= 0   %When on y axis
    theta = pi/2;
elseif x == 0 && y <= 0
    theta = 3*pi/2;
end
end