%% Init
clear all
close all
clc

%% Simulation
% x = [x y z u  v  w  phi  theta  psi  p  q  r]
% u = [T_MR  T_TR  beta_1s  beta_1c]

t_span = [0 10];

x_0 = zeros(12,1);

x_0(5) = 1;

u_0 = [60 -0.15 0 0];

[u_0,fval] = trim(x_0,u_0);


[t,x] = ode45(@(t,x) sys(t,x,u_0),t_span,x_0);


%% Plot

pos_axis    = [0 10 -1e1  1e1];
speed_axis  = [0 10 -5e0  5e0];
rate_axis   = [0 10 -5e0  5e0];
angle_axis  = [0 10 -1e1  1e1];


figure
subplot(2,2,1)
plot(t,x(:,1))
hold on
plot(t,x(:,2))
plot(t,x(:,3))
legend('x','y','z')
title('position')
axis(pos_axis)
grid on

subplot(2,2,2)
plot(t,x(:,4))
hold on
plot(t,x(:,5))
plot(t,x(:,6))
legend('u','v','w')
title('speed')
axis(speed_axis)
grid on

subplot(2,2,3)
plot(t,rad2deg(x(:,7)))
hold on
plot(t,rad2deg(x(:,8)))
plot(t,rad2deg(x(:,9)))
legend('theta','phi','psi')
title('angle')
axis(angle_axis)
grid on

subplot(2,2,4)
plot(t,rad2deg(x(:,10)))
hold on
plot(t,rad2deg(x(:,11)))
plot(t,rad2deg(x(:,12)))
legend('p','q','r')
title('rate')
axis(rate_axis)
grid on