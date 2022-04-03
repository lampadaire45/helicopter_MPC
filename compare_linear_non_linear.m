%% Init
clear all
close all
clc

set(0, 'DefaultLineLineWidth', 1.0);

%% Get dynamics
% u = [theta_0 theta_c]
% x = [x z u w q theta lambda_i]
init_hover_dynamics

%% Simulation conditions
x_0;
u_0;
t_span = [0 10];

% Set step function
u_step_col = @(t,x) u_0 +[0.25*heaviside(t-1)-0.25*heaviside(t-2) 0 ]';
u_step_cyc = @(t,x) u_0 +[0 0.25*heaviside(t-1)-0.25*heaviside(t-2)]';

%% Non-linear simulation
f_step_col = @(t,x) dynamics(t,x,u_step_col,heli_param);
f_step_cyc = @(t,x) dynamics(t,x,u_step_cyc,heli_param);

[t_nl_col,x_nl_col] = ode45(f_step_col,t_span,x_0);
[t_nl_cyc,x_nl_cyc] = ode45(f_step_cyc,t_span,x_0);

%% Plot Non-linear results
% pos_axis     = [0 t_span(2) -5e1   5e1];
% speed_axis   = [0 t_span(2) -1e1   1e1];
% angle_axis   = [0 t_span(2) -5e0   5e0];
% inflow_axis  = [0 t_span(2) -5e-2  5e-2];

pos_axis     = [-inf inf -inf inf];
speed_axis   = [-inf inf -inf inf];
angle_axis   = [-inf inf -inf inf];
inflow_axis  = [-inf inf -inf inf];

axis = [pos_axis;speed_axis;angle_axis;inflow_axis];

plot_single(t_nl_col,x_nl_col,u_step_col,'Non-linear model collective step',axis)
plot_single(t_nl_cyc,x_nl_cyc,u_step_cyc,'Non-linear model cyclic step',axis)

%% Simulate linear system
t_lin = [t_span(1):(t_span(2)-t_span(1))/100:t_span(2)]';
[~,t_lin_col,x_lin_col] = lsim(sysc,fu2u(u_step_col,t_lin,zeros(length(t_lin),length(x_0)))-u_0'.*ones(length(t_lin),size(sysc.B,2)),t_lin);
[~,t_lin_cyc,x_lin_cyc] = lsim(sysc,fu2u(u_step_cyc,t_lin,zeros(length(t_lin),length(x_0)))-u_0'.*ones(length(t_lin),size(sysc.B,2)),t_lin);

%% Plot linear system
plot_single(t_lin_col,x_lin_col+x_0',u_step_col,'Linear model collective step',axis)
plot_single(t_lin_cyc,x_lin_cyc+x_0',u_step_cyc,'Linear model cyclic step',axis)

%% Compare linear and non-linear system
plot_compare_nl2l(t_lin_col,x_lin_col+x_0',u_step_col,t_nl_col,x_nl_col,u_step_col,'Comparison of Linear and Non-Linar System with Collective Step',axis)
plot_compare_nl2l(t_lin_cyc,x_lin_cyc+x_0',u_step_cyc,t_nl_cyc,x_nl_cyc,u_step_cyc,'Comparison of Linear and Non-Linar System with Collective Step',axis)
