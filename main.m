%% Init
clear all
close all
clc

%% Simulation
% x = [x y z u  v  w  phi  theta  psi  p  q  r]
% u = [T_MR  T_TR  beta_1s  beta_1c]

t_span = [0 10];

x_0 = zeros(12,1);

%x_0(5) = 1;

u_0 = [60 -0.15 0 0];

[u_0,fval] = trim(x_0,u_0);

% Set step function
u_step = @(t) ones(size(t))*u_0 +[heaviside(t-5) zeros(size(t)) zeros(size(t)) zeros(size(t))];
f_step = @(t,x) sys(t,x,u_step);

[t,x] = ode45(f_step,t_span,x_0);

%% Plot Non-linear model

pos_axis    = [0 10 -1e1  1e1];
speed_axis  = [0 10 -5e0  5e0];
rate_axis   = [0 10 -5e0  5e0];
angle_axis  = [0 10 -1e1  1e1];

figure
subplot(3,2,1)
plot(t,x(:,1))
hold on
plot(t,x(:,2))
plot(t,x(:,3))
legend('x','y','z')
title('position')
axis(pos_axis)
grid on

subplot(3,2,2)
plot(t,x(:,4))
hold on
plot(t,x(:,5))
plot(t,x(:,6))
legend('u','v','w')
title('speed')
axis(speed_axis)
grid on

subplot(3,2,3)
plot(t,rad2deg(x(:,7)))
hold on
plot(t,rad2deg(x(:,8)))
plot(t,rad2deg(x(:,9)))
legend('theta','phi','psi')
title('angle')
axis(angle_axis)
grid on

subplot(3,2,4)
plot(t,rad2deg(x(:,10)))
hold on
plot(t,rad2deg(x(:,11)))
plot(t,rad2deg(x(:,12)))
legend('p','q','r')
title('rate')
axis(rate_axis)
grid on

subplot(3,1,3)
plot(t,u_step(t)-u_0)
hold on
legend('T_mr','T_tr','\beta_1s','\beta_1c')
title('Command')
grid on

%% Linearize system
[A,B] = lin_sys(x_0,u_0);
sys_lin = ss(A,B,eye(size(A)),zeros(size(B)));

t_lin = [t_span(1):(t_span(2)-t_span(1))/100:t_span(2)]';
[~,t_lin,x_lin] = lsim(sys_lin,u_step(t_lin)-u_0,t_lin);

%% Plot linear system

figure
subplot(3,2,1)
plot(t_lin,x_lin(:,1))
hold on
plot(t_lin,x_lin(:,2))
plot(t_lin,x_lin(:,3))
legend('x','y','z')
title('position')
axis(pos_axis)
grid on

subplot(3,2,2)
plot(t_lin,x_lin(:,4))
hold on
plot(t_lin,x_lin(:,5))
plot(t_lin,x_lin(:,6))
legend('u','v','w')
title('speed')
axis(speed_axis)
grid on

subplot(3,2,3)
plot(t_lin,rad2deg(x_lin(:,7)))
hold on
plot(t_lin,rad2deg(x_lin(:,8)))
plot(t_lin,rad2deg(x_lin(:,9)))
legend('theta','phi','psi')
title('angle')
axis(angle_axis)
grid on

subplot(3,2,4)
plot(t_lin,rad2deg(x_lin(:,10)))
hold on
plot(t_lin,rad2deg(x_lin(:,11)))
plot(t_lin,rad2deg(x_lin(:,12)))
legend('p','q','r')
title('rate')
axis(rate_axis)
grid on

subplot(3,1,3)
plot(t_lin,u_step(t_lin)-u_0)
hold on
legend('T_mr','T_tr','\beta_1s','\beta_1c')
title('Command')
grid on
