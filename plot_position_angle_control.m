function [] = plot_position_angle_control(t,x,input,figure_title,varargin)

% Take care of u if it is a function
u = fu2u(input,t,x);

% Deal with non-exitent axis limits
if ~isempty(varargin)
    ref = varargin{1};
end

figure('name',figure_title)

subplot(3,1,1)
stairs(t,x(:,1))
hold on
if ~isempty(varargin)
    stairs(t,ref(:,1),'--')
end
stairs(t,-x(:,2))
if ~isempty(varargin)
    stairs(t,-ref(:,2),'--')
end
if ~isempty(varargin)
    legend('x','x reference','-z','-z reference')
else
    legend('x','-z')
end
xlabel('Time (s)')
ylabel('Position (m)')
grid on

subplot(3,1,2)
stairs(t,rad2deg(x(:,5)))
hold on
stairs(t,rad2deg(x(:,6)))
legend('q','\theta')
xlabel('Time (s)')
ylabel('Rate and angle (deg/s)')
grid on

subplot(3,1,3)
stairs(t,u)
legend('Collective','Cyclic')
xlabel('Time (s)')
ylabel('Command (% of maximum)')
grid on
end