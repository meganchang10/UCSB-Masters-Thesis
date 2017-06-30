close all
clear all
clc

m = 200;
iterations = 40;
D = 0.5;
phi_eqn = @(x,y) max(abs(x),abs(y)) - 0.1*sqrt(2);    %Square
d0 = 0.005;


tfinal = 0.0252;
epsilon = 0.7;
[times_through_loopA,phi,FramesA] = Snowflake_movie(m,iterations,D,d0,epsilon,phi_eqn,tfinal);

tfinal = 0.025;
epsilon = 0.6;
[times_through_loopB,phi,FramesB] = Snowflake_movie(m,iterations,D,d0,epsilon,phi_eqn,tfinal);

tfinal = 0.0502;
epsilon = 0.3;
[times_through_loopC,phi,FramesC] = Snowflake_movie(m,iterations,D,d0,epsilon,phi_eqn,tfinal);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% video = VideoWriter('Snowflake_A.avi','Uncompressed AVI');
% open(video)
% 
% for i = 1:155
% writeVideo(video,FramesA(i));
% end
% for i = 190:296
% writeVideo(video,FramesA(i));
% end
% for i = 325:449
% writeVideo(video,FramesA(i));
% end
% 
% close(video)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% video = VideoWriter('Snowflake_B.avi','Uncompressed AVI');
% open(video)
% 
% for i = 1:330
% writeVideo(video,FramesB(i));
% end
% 
% close(video)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% video = VideoWriter('Snowflake_C.avi','Uncompressed AVI');
% open(video)
% 
% for i = 1:135
% writeVideo(video,FramesC(i))
% end
% for i = 185:220
% writeVideo(video,FramesC(i))
% end
% for i = 235:300
% writeVideo(video,FramesC(i))
% end
% for i = 320:450
% writeVideo(video,FramesC(i))
% end
% 
% close(video)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
