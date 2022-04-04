function [] = plot_compare_nl2l(t_lin,x_lin,input_lin,t_nl,x_nl,input_nl,figure_title,varargin)

% Take care of u if it is a function
u = fu2u(input_nl,t_nl,x_nl);

% Take care of u if it is a function
u_lin = fu2u(input_lin,t_lin,x_lin);

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
stairs(t_nl,x_nl(:,1))
hold on
stairs(t_nl,-x_nl(:,2))
stairs(t_lin,x_lin(:,1),'--')
stairs(t_lin,-x_lin(:,2),'--')
legend('x non-linear','-z non-linear','x linear','-z linear')
xlabel('Time (s)')
ylabel('Position (m)')
grid on
axis(pos_axis)

subplot(3,2,2)
stairs(t_nl,x_nl(:,3))
hold on
stairs(t_nl,-x_nl(:,4))
stairs(t_lin,x_lin(:,3),'--')
stairs(t_lin,-x_lin(:,4),'--')
legend('u non-linear','-w non-linear','u linear','-w linear')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
grid on
axis(speed_axis)

subplot(3,2,3)
stairs(t_nl,rad2deg(x_nl(:,5)))
hold on
stairs(t_nl,rad2deg(x_nl(:,6)))
stairs(t_lin,rad2deg(x_lin(:,5)),'--')
stairs(t_lin,rad2deg(x_lin(:,6)),'--')
legend('q non-linear','\theta non-linear','q linear','\theta linear')
xlabel('Time (s)')
ylabel('Rate and angle (deg/s)')
grid on
axis(angle_axis)

subplot(3,2,4)
stairs(t_nl,x_nl(:,7))
hold on
stairs(t_lin,x_lin(:,7),'--')
legend('Inflow non-linear','Inflow linear')
xlabel('Time (s)')
ylabel('Inflow')
grid on
axis(inflow_axis)

subplot(3,1,3)
stairs(t_nl,u)
hold on
stairs(t_lin,u_lin,'--')
legend('Collective non-linear','Cyclic non-linear','Collective linear','Cyclic linear')
xlabel('Time (s)')
ylabel('Command (% of maximum)')
grid on
end