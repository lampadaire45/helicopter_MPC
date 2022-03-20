function [] = plot_compare_nl2l(t_lin,x_lin,input_lin,t_nl,x_nl,input_nl,figure_title,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Take care of u if it is a function
u_lin = fu2u(input_lin,t_lin,x_lin);

u_nl = fu2u(input_nl,t_nl,x_nl);

% Deal with non-exitent axis limits
if isempty(varargin)
    pos_axis = [-inf inf -inf inf];
    speed_axis = [-inf inf -inf inf];
    rate_axis = [-inf inf -inf inf];
    angle_axis = [-inf inf -inf inf];
else
    axis_array = varargin{1};
    pos_axis = axis_array(1,:);
    speed_axis = axis_array(2,:);
    rate_axis = axis_array(3,:);
    angle_axis = axis_array(4,:);
end

% Separate x vector
x_pos_lin = x_lin(:,[1:3]);
x_speed_lin = x_lin(:,[4:6]);
x_angle_lin = rad2deg(x_lin(:,[7:9]));
x_rate_lin = rad2deg(x_lin(:,[10:12]));

x_pos_nl = x_nl(:,[1:3]);
x_speed_nl = x_nl(:,[4:6]);
x_angle_nl = rad2deg(x_nl(:,[7:9]));
x_rate_nl = rad2deg(x_nl(:,[10:12]));

% Plot
figure('name',figure_title)
subplot(3,2,1)
plot(t_lin,x_pos_lin)
hold on
plot(t_nl,x_pos_nl)
legend('x (linear)','y (linear)','z (linear)','x (non-linear)','y (non-linear)','z (non-linear)')
%title('Position')
xlabel('Time (s)');ylabel('Position (m)')
axis(pos_axis)
grid on

subplot(3,2,2)
plot(t_lin,x_speed_lin)
hold on
plot(t_nl,x_speed_nl)
legend('u (linear)','v (linear)','w (linear)','u (non-linear)','v (non-linear)','w (non-linear)')
%title('speed')
xlabel('Time (s)');ylabel('Speed (m/s)')
axis(speed_axis)
grid on

subplot(3,2,3)
plot(t_lin,x_angle_lin)
hold on
plot(t_nl,x_angle_nl)
legend('\theta (linear)','\phi (linear)','\psi (linear)','\theta (non-linear)','\phi (non-linear)','\psi (non-linear)')
xlabel('Time (s)');ylabel('Attitude (deg)')
%title('angle')
axis(angle_axis)
grid on

subplot(3,2,4)
plot(t_lin,x_rate_lin)
hold on
plot(t_nl,x_rate_nl)
legend('p (linear)','q (linear)','r (linear)','p (non-linear)','q (non-linear)','r (non-linear)')
xlabel('Time (s)');ylabel('Rate (deg/s)')
%title('rate')
axis(rate_axis)
grid on

subplot(3,1,3)
yyaxis left
plot(t_lin,u_lin(:,[1,2]))
hold on
plot(t_nl,u_nl(:,[1,2]))
ylabel('Thrust command (N)')

yyaxis right
plot(t_lin,rad2deg(u_lin(:,[3,4])))
hold on
plot(t_nl,rad2deg(u_nl(:,[3,4])))
ylabel('Angle command (deg)')
legend('T_{mr} (linear)','T_{mr} (non-linear)','T_{tr} (linear)', ...
    'T_{tr} (non-linear)','\beta_{1s} (linear)','\beta_{1s} (non-linear)', ...
    '\beta_{1c} (linear)','\beta_{1c} (non-linear)')
xlabel('Time (s)');
%title('Command')
grid on


end