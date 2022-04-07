%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare linear to non-linear model for longitudinal
% helicopter model
%
% By: Frédéric Larocque
%
% Last modified: 6/04/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Init
clear all
close all
clc
addpath('functions')

set(0, 'DefaultLineLineWidth', 1.0);

%% Get dynamics
% u = [theta_0 theta_c]
% x = [x z u w q theta lambda_i]
init_hover_dynamics

fprintf('Initialized linear and non-linear system dynamics\n')
%% Simulation conditions
x_0 = x_0_eq;
u_0 = u_0_eq;
t_span = [0 10];

% Set step function
u_step_col = @(t,x) u_0 +[0.25*heaviside(t-1)-0.25*heaviside(t-2) 0 ]';
u_step_cyc = @(t,x) u_0 +[0 0.25*heaviside(t-1)-0.25*heaviside(t-2)]';

%% Non-linear simulation
% Set dynamics of system + input
f_step_col = @(t,x) dynamics(t,x,u_step_col,heli_param);
f_step_cyc = @(t,x) dynamics(t,x,u_step_cyc,heli_param);

% Solve using ode45
[t_nl_col,x_nl_col] = ode45(f_step_col,t_span,x_0);
[t_nl_cyc,x_nl_cyc] = ode45(f_step_cyc,t_span,x_0);

fprintf('Simulated non-linear system\n')
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

fprintf('Simulated linear system\n')
%% Plot linear system
plot_single(t_lin_col,x_lin_col+x_0',u_step_col,'Linear model collective step',axis)
plot_single(t_lin_cyc,x_lin_cyc+x_0',u_step_cyc,'Linear model cyclic step',axis)

%% Compare linear and non-linear system
plot_compare_nl2l(t_lin_col,x_lin_col+x_0',u_step_col,t_nl_col,x_nl_col,u_step_col,'Comparison of Linear and Non-Linar System with Collective Step',axis)
plot_compare_nl2l(t_lin_cyc,x_lin_cyc+x_0',u_step_cyc,t_nl_cyc,x_nl_cyc,u_step_cyc,'Comparison of Linear and Non-Linar System with Cyclic Step',axis)

%% Plot for report

figure ('name','Comparison of linear and non-linear system for velocity and rate, cyclic step')
subplot(2,1,1)
stairs(t_nl_cyc,x_nl_cyc(:,3))
hold on
stairs(t_nl_cyc,-x_nl_cyc(:,4))
stairs(t_lin_cyc,x_lin_cyc(:,3))
stairs(t_lin_cyc,-x_lin_cyc(:,4))
legend('u non-linear','-w non-linear','u linear','-w linear')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
grid on

subplot(2,1,2)
stairs(t_nl_cyc,rad2deg(x_nl_cyc(:,5)))
hold on
stairs(t_nl_cyc,rad2deg(x_nl_cyc(:,6)))
stairs(t_lin_cyc,rad2deg(x_lin_cyc(:,5)))
stairs(t_lin_cyc,rad2deg(x_lin_cyc(:,6)))
legend('q non-linear','\theta non-linear','q linear','\theta linear')
xlabel('Time (s)')
ylabel('Rate and angle (deg/s)')
grid on
