%% Init
clear all
close all
clc

%% Simulation

t_span = [0 10];

x_0 = zeros(9,1);

%u_0 = [68.67 -0.1582 0 0];
u_0 = [72 -0.172 0 0];

[t,x] = ode45(@(t,x) sys(t,x,u_0),t_span,x_0);


%% Plot
figure
subplot(2,2,1)
plot(t,x(:,1))
hold on
plot(t,x(:,2))
plot(t,x(:,3))
legend('u','v','w')
title('position')

subplot(2,2,2)
plot(t,x(:,4))
hold on
plot(t,x(:,5))
plot(t,x(:,6))
legend('theta','phi','psi')
title('angle')

subplot(2,1,2)
plot(t,x(:,7))
hold on
plot(t,x(:,8))
plot(t,x(:,9))
legend('p','q','r')
title('rate')