%% Init
clear all
close all
clc

set(0, 'DefaultLineLineWidth', 1.0);
%% Simulation Conditions
% x = [x y z u  v  w  phi  theta  psi  p  q  r]
% u = [T_MR  T_TR  beta_1s  beta_1c]

t_span = [0 10];
x_0 = zeros(12,1);
u_0 = [60 -0.15 0 0]';
[u_0,fval] = trim(x_0,u_0);

% Set step function
u_step_col = @(t,x) u_0 +[10*heaviside(t-1)-10*heaviside(t-2) 0 0 0]';
u_step_cyc = @(t,x) u_0 +[0 0 deg2rad(1*heaviside(t-1)-1*heaviside(t-2)) 0]';

%% Non-linear simulation
f_step_col = @(t,x) sys(t,x,u_step_col);
f_step_cyc = @(t,x) sys(t,x,u_step_cyc);

[t_nl_col,x_nl_col] = ode45(f_step_col,t_span,x_0);
[t_nl_cyc,x_nl_cyc] = ode45(f_step_cyc,t_span,x_0);

%% Plot Non-linear results
pos_axis    = [0 10 -1e1  1e1];
speed_axis  = [0 10 -5e0  5e0];
rate_axis   = [0 10 -5e0  5e0];
angle_axis  = [0 10 -1e1  1e1];
axis = [pos_axis;speed_axis;rate_axis;angle_axis];

plot_single(t_nl_col,x_nl_col,u_step_col,'Non-linear model collective step')
plot_single(t_nl_cyc,x_nl_cyc,u_step_cyc,'Non-linear model cyclic step')
% Highly coupled system

visualize_helicopter_trajectory_rotating(x_nl_col,[],0.01)

%% Linearize system in hover
% x = [x y z u  v  w  phi  theta  psi  p  q  r]
% u = [T_MR  T_TR  beta_1s  beta_1c]
x_0 = zeros(12,1);
u_0 = [60 -0.15 0 0]';
[u_0,fval] = trim(x_0,u_0);

[A,B] = lin_sys(x_0,u_0);
sys_lin = ss(A,B,eye(size(A)),zeros(size(B)));
t_lin = [t_span(1):(t_span(2)-t_span(1))/100:t_span(2)]';

%% Simulate linear
[~,t_lin_col,x_lin_col] = lsim(sys_lin,fu2u(u_step_col,t_lin,zeros(length(t_lin),length(x_0)))-u_0,t_lin);
[~,t_lin_cyc,x_lin_cyc] = lsim(sys_lin,fu2u(u_step_cyc,t_lin,zeros(length(t_lin),length(x_0)))-u_0,t_lin);

%% Plot linear system
plot_single(t_lin_col,x_lin_col,u_step_col,'Linear model collective step')
plot_single(t_lin_cyc,x_lin_cyc,u_step_cyc,'Linear model cyclic step')

%% Compare linear and non-linear system
plot_compare_nl2l(t_lin_col,x_lin_col,u_step_col,t_nl_col,x_nl_col,u_step_col,'Comparison of Linear and Non-Linar System with Collective Step',axis)
plot_compare_nl2l(t_lin_cyc,x_lin_cyc,u_step_cyc,t_nl_cyc,x_nl_cyc,u_step_cyc,'Comparison of Linear and Non-Linar System with Collective Step',axis)

%% Check controlability of system

if rank(ctrb(sys_lin.A,sys_lin.B)) == size(sys_lin.A,1)
    fprintf("System is controllable")
else
    fprintf("System not controllable")
end

%% Discretize system
dt = 0.05;

sysd = c2d(sys_lin,dt);

%% LQR controller
% x = [x y z u  v  w  phi  theta  psi  p  q  r]
x_0 = [5 5 0 0  0  0  0    0      0    0  0  0]';

% u = [T_MR  T_TR  beta_1s  beta_1c]
u_0 = [0     0     0        0]';

t_span = [0 15];
t = [t_span(1):dt:t_span(2)]';

Q = 1E3*diag([1 1 1 1 1 1 1 1 1 1 1 1]);
R = diag([1 1 1 1]);

[K,S,e] = dlqr(sysd.A,sysd.B,Q,R,[]);

sysd_LQR_cl = ss(sysd.A-sysd.B*K,sysd.B,sysd.C,sysd.D,dt);
[~,t_LQR,x_LQR] = lsim(sysd_LQR_cl,zeros(length(t),size(sysd.B,2)),t,x_0);

u_LQR = zeros(length(t_LQR),size(sysd.B,2));
for i=1:length(t_LQR)
    u_LQR(i,:) = K*x_LQR(i,:)';
end

plot_single(t_LQR,x_LQR,u_LQR,'LQR')

visualize_helicopter_trajectory_rotating(x_LQR,[],0.01)