function [] = plot_single(t,x,input,figure_title,varargin)

% Take care of u if it is a function
u = fu2u(input,t,x);

% Deal with non-exitent axis limits
if isempty(varargin)
    pos_axis = [-inf inf -inf inf];
    speed_axis = [-inf inf -inf inf];
    angle_axis = [-inf inf -inf inf];
    inflow_axis = [-inf inf -inf inf];
else
    pos_axis = varargin{1}(1,:);
    speed_axis = varargin{1}(2,:);
    angle_axis = varargin{1}(3,:);
    inflow_axis = varargin{1}(4,:);
end

figure('name',figure_title)

subplot(3,2,1)
stairs(t,x(:,1))
hold on
stairs(t,-x(:,2))
legend('x','-z')
xlabel('Time (s)')
ylabel('Position (m)')
grid on
axis(pos_axis)

subplot(3,2,2)
stairs(t,x(:,3))
hold on
stairs(t,-x(:,4))
legend('u','-w')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
grid on
axis(speed_axis)

subplot(3,2,3)
stairs(t,rad2deg(x(:,5)))
hold on
stairs(t,rad2deg(x(:,6)))
legend('q','\theta')
xlabel('Time (s)')
ylabel('Rate and angle (deg/s)')
grid on
axis(angle_axis)

subplot(3,2,4)
stairs(t,x(:,7))
legend('Inflow')
xlabel('Time (s)')
ylabel('Inflow')
grid on
axis(inflow_axis)

subplot(3,1,3)
stairs(t,u)
legend('Collective','Cyclic')
xlabel('Time (s)')
ylabel('Command (% of maximum)')
grid on
end