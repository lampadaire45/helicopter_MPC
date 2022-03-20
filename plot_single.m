function [] = plot_single(t,x,input,figure_title,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Take care of u if it is a function
u = fu2u(input,t,x);

% Deal with non-exitent axis limits
if isempty(varargin)
    pos_axis = [-inf inf -inf inf];
    speed_axis = [-inf inf -inf inf];
    rate_axis = [-inf inf -inf inf];
    angle_axis = [-inf inf -inf inf];
else
    pos_axis = varargin(1,:);
    speed_axis = varargin(2,:);
    rate_axis = varargin(3,:);
    angle_axis = varargin(4,:);
end

% Separate x vector
x_pos = x(:,[1:3]);
x_speed = x(:,[4:6]);
x_angle = rad2deg(x(:,[7:9]));
x_rate = rad2deg(x(:,[10:12]));

% Plot
figure('name',figure_title)
subplot(3,2,1)
plot(t,x_pos)
hold on
legend('x','y','z')
%title('Position')
xlabel('Time (s)');ylabel('Position (m)')
axis(pos_axis)
grid on

subplot(3,2,2)
plot(t,x_speed)
hold on
legend('u','v','w')
%title('speed')
xlabel('Time (s)');ylabel('Speed (m/s)')
axis(speed_axis)
grid on

subplot(3,2,3)
plot(t,x_angle)
hold on
legend('\theta','\phi','\psi')
xlabel('Time (s)');ylabel('Attitude (deg)')
%title('angle')
axis(angle_axis)
grid on

subplot(3,2,4)
plot(t,x_rate)
hold on
legend('p','q','r')
xlabel('Time (s)');ylabel('Rate (deg/s)')
%title('rate')
axis(rate_axis)
grid on

subplot(3,1,3)
yyaxis left
plot(t,u(:,[1,2]))
ylabel('Thrust command (N)')
hold on
yyaxis right
plot(t,rad2deg(u(:,[3,4])))
ylabel('Angle command (deg)')
legend('T_{mr}','T_{tr}','\beta_{1s}','\beta_{1c}')
xlabel('Time (s)');
%title('Command')
grid on


end