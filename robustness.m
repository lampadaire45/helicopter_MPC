%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robustness test on MPC on helicopter model
%
% Script to analyse effect of change of model on MPC controller
%
% By: Frédéric Larocque
%
% Last modified: 7/04/2022
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
sysd_real = sysd;
init_hover_dynamics_robustness_test;

fprintf('Initialized real and modified system dynamics\n')

%% LTI system definition (modified system)
% x = [x    z   u   w   q   theta lambda_i]'
x_0 = [0.1  0   0   0   0   0     0     ]'; % Initial condition inside Xf

% u = [theta_0 theta_c]'
u_0 = [0       0]';

% d = [delta_u delta_w]'
% y = [x y q w]'

%Definition of the LTI system
LTI.A=sysd_real.A; 
LTI.B=sysd_real.B;
LTI.C=[eye(size(LTI.A))];
LTI.x0=x_0;

%Definition of system dimension
dim.nx=size(LTI.A,1);     %state  dimension
dim.nu=size(LTI.B,2);     %input  dimension
dim.ny=size(LTI.C,1);     %output dimension
dim.N=25;      %horizon

%Definition of quadratic cost function
%           x = [    x    z    u     w    q    theta   lambda_i]
weight.Q = 1E1*diag([1E0  1E0  1E0   1E0  1E0  1E0     1E-4]);
%           u = [theta_0 theta_c]
weight.R = diag([1       1]);
weight.beta = 1;

[K,weight.P,e] = dlqr(LTI.A,LTI.B,weight.Q,weight.R,[]);
K = -K;

% Constraints
%   u = [theta_0 theta_c] % percentage of max value
u_lim = [0.25    1]';
%   x = [x    z    u   w   q            theta          lambda_i]
x_lim = [1000 1000 5   5   deg2rad(5)   deg2rad(15)    1000]';
x_lim_vec = repmat(x_lim,[dim.N+1,1]);

fprintf('General MPC Parameters initialized\n')

%% MPC regulation for modified system
% x = [x    z   u   w   q   theta lambda_i]'
x_0 = [1  0   0   0   0   0     0     ]'; % Initial condition inside Xf
LTI.x0 = x_0;

% Generation of prediction model 
predmod=predmodgen(LTI,dim);            
[H,h]=costgen(predmod,weight,dim);

% Time vectors
t_span = [0 10];
t = [t_span(1):sysd.Ts:t_span(2)]';

T=length(t);                  %simulation horizon
x=zeros(dim.nx,T+1);
u_rec=zeros(dim.nu,T);
x(:,1)=LTI.x0;

fprintf('Starting modified MPC simulation\n')
% Receding horizon implementation
for k=1:T
    fprintf('Time is %2.2f \n',t(k))

    x_0_MPC=x(:,k);

    % Solve problem with CVX
    cvx_begin quiet
        variable uN(dim.nu*dim.N)
        minimize(0.5*uN'*H*uN+(h*x_0_MPC)'*uN)
        
        % input constraints
        uN <=  repmat(u_lim,[dim.N,1]);
        uN >= -repmat(u_lim,[dim.N,1]);
        % state constraints
        predmod.S*uN <= -predmod.T*x_0_MPC + x_lim_vec;
        predmod.S*uN >= -predmod.T*x_0_MPC - x_lim_vec; 

    cvx_end
    % Select the first input only
    u_rec(:,k)=uN(1:1*dim.nu);

    % Compute the state/output evolution
    x(:,k+1)=sysd.A*x_0_MPC + sysd.B*u_rec(:,k); %modified system
    
end
u_robust = u_rec;
x_robust = x(:,1:end-1);
t_robust = t;
clc
fprintf('Simulation finished\n')

plot_position_angle_control(t_robust,x_robust',u_robust','Robust')


%% Real system
% x = [x    z   u   w   q   theta lambda_i]'
x_0 = [1  0   0   0   0   0     0     ]'; % Initial condition inside Xf
LTI.x0 = x_0;
x(:,1)=LTI.x0;

fprintf('Starting normal MPC simulation\n')
% Receding horizon implementation
for k=1:T
    fprintf('Time is %2.2f \n',t(k))

    x_0_MPC=x(:,k);

    % Solve problem with CVX
    cvx_begin quiet
        variable uN(dim.nu*dim.N)
        minimize(0.5*uN'*H*uN+(h*x_0_MPC)'*uN)
        
        % input constraints
        uN <=  repmat(u_lim,[dim.N,1]);
        uN >= -repmat(u_lim,[dim.N,1]);
        % state constraints
        predmod.S*uN <= -predmod.T*x_0_MPC + x_lim_vec;
        predmod.S*uN >= -predmod.T*x_0_MPC - x_lim_vec; 

    cvx_end
    % Select the first input only
    u_rec(:,k)=uN(1:1*dim.nu);

    % Compute the state/output evolution
    x(:,k+1)=LTI.A*x_0_MPC + LTI.B*u_rec(:,k); %modified system
    
end
u_norm = u_rec;
x_norm = x(:,1:end-1);
t_norm = t;
clc
fprintf('Simulation finished\n')

plot_position_angle_control(t_norm,x_norm',u_norm','Robust')

%% Comparison
figure('name','Comparison of normal and modified system')
subplot(2,1,1)
stairs(t_norm,x_norm(1,:))
hold on
stairs(t_robust,x_robust(1,:),'--')
legend('Normal','Modified System')
xlabel('Time (s)')
ylabel('Position (m)')
grid on

subplot(2,1,2)
stairs(t_norm,u_norm(2,:))
hold on
stairs(t_robust,u_robust(2,:),'--')
legend('Normal','Modified System')
xlabel('Time (s)')
ylabel('Cyclic Command (% of maximum)')
grid on
