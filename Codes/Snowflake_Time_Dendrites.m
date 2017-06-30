% close all
% clear all
% clc


m = 200;
iterations = 40;

phi_eqn = @(x,y) max(abs(x),abs(y)) - 0.1*sqrt(2);    %Square

time_steps_or_nah = 1; %if you want figure of zero level set every couple seconds
% time_steps_or_nah = 0; %if you only need the one figure with each iteration on the same plot

d0 = 0.005;


dt_plot = 0.005;
tfinal = 0.0252;

figure_number = 1;
epsilon = 0.7;
[times_through_loopA,phiA] = Snowflake_fun(m,iterations,D,d0,epsilon,dt_plot,phi_eqn,tfinal,figure_number,time_steps_or_nah);

figure_number = 7;
epsilon = 0.6;
[times_through_loopB,phiB] = Snowflake_fun(m,iterations,D,d0,epsilon,dt_plot,phi_eqn,tfinal,figure_number,time_steps_or_nah);


dt_plot = 0.010;
tfinal = 0.0502615;

figure_number = 13;
epsilon = 0.3;
[times_through_loopC, phiC] = Snowflake_fun(m,iterations,D,d0,epsilon,dt_plot,phi_eqn,tfinal,figure_number,time_steps_or_nah);
