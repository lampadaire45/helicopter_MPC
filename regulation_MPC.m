%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regulation MPC on helicopter model
%
% Script to analyse effect of N, effect of Q and R values and comparison to
% LQR controller
%
% By: Frédéric Larocque
%
% Last modified: 6/04/2022
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

fprintf('Initialized linear and non-linear system dynamics\n')

%% LTI system definition
% x = [x    z   u   w   q   theta lambda_i]'
x_0 = [0.1  0   0   0   0   0     0     ]'; % Initial condition inside Xf

% u = [theta_0 theta_c]'
u_0 = [0       0]';

% d = [delta_u delta_w]'
% y = [x y q w]'

%Definition of the LTI system
LTI.A=sysd.A; 
LTI.B=sysd.B;
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
%% MPC regulation in Xf
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

fprintf('Starting regulation MPC simulation\n')
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
    x(:,k+1)=LTI.A*x_0_MPC + LTI.B*u_rec(:,k);
    
end
u_MPC_inXf = u_rec;
x_MPC_inXf = x(:,1:end-1);
t_MPC_inXf = t;
clc
fprintf('Simulation finished\n')

% Plots
plot_position_control(t_MPC_inXf,x_MPC_inXf',u_MPC_inXf','Regulation MPC in Xf')

%% LQR simulation in Xf
t_LQR_inXf = t;
x_LQR_inXf = zeros(length(x_0),length(t_LQR_inXf));
x_LQR_inXf(:,1) = x_0;
u_LQR_inXf = zeros(length(u_0),length(t_LQR_inXf));

fprintf('Starting regulation LQR simulation\n')
for i=1:length(t)-1
    fprintf('Time is %2.2f \n',t(i))
    % Calculate u
    u_LQR_inXf(:,i) = K*x_LQR_inXf(:,i);

    x_LQR_inXf(:,i+1) = LTI.A*x_LQR_inXf(:,i)+LTI.B*u_LQR_inXf(:,i);

end
clc
fprintf('Simulation finished\n')
% Plots
plot_position_control(t_LQR_inXf,x_LQR_inXf',u_LQR_inXf','Regulation LQR in Xf')

%% Comparison of LQR and MPC in Xf
figure('name','Comparison of LQR and MPC Regulation when in Xf')
subplot(2,1,1)
stairs(t_MPC_inXf,x_MPC_inXf(1,:))
hold on
stairs(t_LQR_inXf,x_LQR_inXf(1,:),'--')
legend('MPC','LQR')
xlabel('Time (s)')
ylabel('Position (m)')
grid on

subplot(2,1,2)
stairs(t_MPC_inXf,u_MPC_inXf(2,:))
hold on
stairs(t_LQR_inXf,u_LQR_inXf(2,:),'--')
legend('MPC','LQR')
xlabel('Time (s)')
ylabel('Cyclic Command (% of maximum)')
grid on

%% MPC Outside Xf
% x = [x    z   u   w   q   theta lambda_i]'
x_0 = [2  0   0   0   0   0     0     ]'; % Initial condition outside Xf
LTI.x0 = x_0;
x(:,1)=LTI.x0;

fprintf('Starting regulation MPC simulation\n')
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
    x(:,k+1)=LTI.A*x_0_MPC + LTI.B*u_rec(:,k);
    
end
u_MPC_outXf = u_rec;
x_MPC_outXf = x(:,1:end-1);
t_MPC_outXf = t;
clc
fprintf('Simulation finished\n')

% Plots
plot_position_control(t_MPC_outXf,x_MPC_outXf',u_MPC_outXf','Regulation MPC out of Xf')
%% LQR Outside Xf
t_LQR_outXf = t;
x_LQR_outXf = zeros(length(x_0),length(t_LQR_outXf));
x_LQR_outXf(:,1) = x_0;
u_LQR_outXf = zeros(length(u_0),length(t_LQR_outXf));

fprintf('Starting regulation LQR simulation\n')
for i=1:length(t)-1
    fprintf('Time is %2.2f \n',t(i))
    % Calculate u
    u_LQR_outXf(:,i) = K*x_LQR_outXf(:,i);

    x_LQR_outXf(:,i+1) = LTI.A*x_LQR_outXf(:,i)+LTI.B*u_LQR_outXf(:,i);

end
clc
fprintf('Simulation finished\n')
% Plots
plot_position_control(t_LQR_outXf,x_LQR_outXf',u_LQR_outXf','Regulation LQR out of Xf')

%% Comparison of LQR and MPC out of Xf
figure('name','Comparison of LQR and MPC Regulation when out of Xf')
subplot(2,1,1)
stairs(t_MPC_outXf,x_MPC_outXf(1,:))
hold on
stairs(t_LQR_outXf,x_LQR_outXf(1,:),'--')
legend('MPC','LQR')
xlabel('Time (s)')
ylabel('Position (m)')
grid on

subplot(2,1,2)
stairs(t_MPC_outXf,u_MPC_outXf(2,:))
hold on
stairs(t_LQR_outXf,u_LQR_outXf(2,:),'--')
legend('MPC','LQR')
xlabel('Time (s)')
ylabel('Cyclic Command (% of maximum)')
grid on

%% MPC Changing initial condition
x_list = [0.1 1 10 45];
x_results = cell(length(x_list),3);

% Time vectors
t_span = [0 20];
t = [t_span(1):sysd.Ts:t_span(2)]';

T=length(t);                  %simulation horizon
x=zeros(dim.nx,T+1);
u_rec=zeros(dim.nu,T);
x(:,1)=LTI.x0;


LTI.x0 = x_0;

for i=1:length(x_list)
    % x = [x    z   u   w   q   theta lambda_i]'
    x_0 = [x_list(i)  0   0   0   0   0     0     ]';
    LTI.x0 = x_0;
    x(:,1)=LTI.x0;

    fprintf('Starting regulation MPC simulation x:%2.0f\n',LTI.x0(1))
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
        x(:,k+1)=LTI.A*x_0_MPC + LTI.B*u_rec(:,k);
        
    end
    x_results{i,1} = u_rec;
    x_results{i,2} = x(:,1:end-1);
    x_results{i,3} = t;
    clc
    fprintf('Simulation finished\n')


end

figure('name','Comparison of value of x for MPC regulation')
for i=1:length(x_list)
    subplot(2,1,1)
    stairs(x_results{i,3},x_results{i,2}(1,:))
    hold on
    subplot(2,1,2)
    stairs(x_results{i,3},x_results{i,1}(2,:))
    hold on
end

subplot(2,1,1)
xlabel('Time (s)')
ylabel('Position (m)')
grid on
legend(cellstr(num2str(x_list', 'x = %-d')))
subplot(2,1,2)
xlabel('Time (s)')
ylabel('Cyclic Command (% of maximum)')
grid on
legend(cellstr(num2str(x_list', 'x = %-d')))

%% Changing N
N_list = [5 10 15 25];

N_results = cell(length(N_list),3);

% x = [x    z   u   w   q   theta lambda_i]'
x_0 = [2  0   0   0   0   0     0     ]'; % Initial condition inside Xf
LTI.x0 = x_0;

for i=1:length(N_list)
    % Generation of prediction model 
    dim.N = N_list(i);
    predmod=predmodgen(LTI,dim);            
    [H,h]=costgen(predmod,weight,dim);
    x_lim_vec = repmat(x_lim,[dim.N+1,1]);
    x(:,1)=LTI.x0;

    fprintf('Starting regulation MPC simulation N:%2.0f\n',dim.N)
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
        x(:,k+1)=LTI.A*x_0_MPC + LTI.B*u_rec(:,k);
        
    end
    N_results{i,1} = u_rec;
    N_results{i,2} = x(:,1:end-1);
    N_results{i,3} = t;
    clc
    fprintf('Simulation finished\n')


end

figure('name','Comparison of value of N for MPC regulation')
for i=1:length(N_list)
    subplot(2,1,1)
    stairs(N_results{i,3},N_results{i,2}(1,:))
    hold on
    subplot(2,1,2)
    stairs(N_results{i,3},N_results{i,1}(2,:))
    hold on
end

subplot(2,1,1)
xlabel('Time (s)')
ylabel('Position (m)')
grid on
legend(cellstr(num2str(N_list', 'N = %-d')))
subplot(2,1,2)
xlabel('Time (s)')
ylabel('Cyclic Command (% of maximum)')
grid on
legend(cellstr(num2str(N_list', 'N = %-d')))

%% Changing Q
% x = [x    z   u   w   q   theta lambda_i]'
x_0 = [0    1.5   0   0   0   0     0     ]'; % Initial condition inside Xf
LTI.x0 = x_0;

% Generation of prediction model 
dim.N = 25;
predmod=predmodgen(LTI,dim);            
[H,h]=costgen(predmod,weight,dim);
x_lim_vec = repmat(x_lim,[dim.N+1,1]);

%Definition of quadratic cost function
%           x = [    x    z    u     w    q    theta   lambda_i]
weight.Q = 1E1*diag([1E0  1E0  1E0   1E0  1E0  1E0     1E-4]);
%           u = [theta_0 theta_c]
weight.R = diag([1       1]);
weight.beta = 1;

[K,weight.P,e] = dlqr(LTI.A,LTI.B,weight.Q,weight.R,[]);
K = -K;

Q_list = [1E-4 1E-2 1E1];
Q_results = cell(length(Q_list),3);

for i=1:length(Q_list)
    % Generation of cost model 
    weight.Q = Q_list(i)*diag([1E0  1E0  1E0   1E0  1E0  1E0     1E-4]);
    [K,weight.P,e] = dlqr(LTI.A,LTI.B,weight.Q,weight.R,[]);
    K = -K;
    [H,h]=costgen(predmod,weight,dim);
    x(:,1)=LTI.x0;

    fprintf('Starting regulation MPC simulation Q:%2.1f\n',Q_list(i))
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
        x(:,k+1)=LTI.A*x_0_MPC + LTI.B*u_rec(:,k);
        
    end
    Q_results{i,1} = u_rec;
    Q_results{i,2} = x(:,1:end-1);
    Q_results{i,3} = t;
    clc
    fprintf('Simulation finished\n')


end

figure('name','Comparison of value of Q for MPC regulation')
for i=1:length(Q_list)
    subplot(2,1,1)
    stairs(Q_results{i,3},Q_results{i,2}(2,:))
    hold on
    subplot(2,1,2)
    stairs(Q_results{i,3},Q_results{i,1}(1,:))
    hold on
end

subplot(2,1,1)
xlabel('Time (s)')
ylabel('Position (m)')
grid on
legend(cellstr(num2str(Q_list', 'Q = %-1.0e')))
subplot(2,1,2)
xlabel('Time (s)')
ylabel('Collective Command (% of maximum)')
grid on
legend(cellstr(num2str(Q_list', 'Q = %-1.0e')))
%% Changing R

% x = [x    z   u   w   q   theta lambda_i]'
x_0 = [0    1.5   0   0   0   0     0     ]'; % Initial condition inside Xf
LTI.x0 = x_0;

% Generation of prediction model 
dim.N = 25;
predmod=predmodgen(LTI,dim);            
[H,h]=costgen(predmod,weight,dim);
x_lim_vec = repmat(x_lim,[dim.N+1,1]);

%Definition of quadratic cost function
%           x = [    x    z    u     w    q    theta   lambda_i]
weight.Q = 1E1*diag([1E0  1E0  1E0   1E0  1E0  1E0     1E-4]);
%           u = [theta_0 theta_c]
weight.R = diag([1       1]);
weight.beta = 1;

[K,weight.P,e] = dlqr(LTI.A,LTI.B,weight.Q,weight.R,[]);
K = -K;

R_list = [1E-2 1E2 1E4];
R_results = cell(length(R_list),3);

for i=1:length(R_list)
    % Generation of cost model 
    weight.R = R_list(i)*diag([1       1]);
    [K,weight.P,e] = dlqr(LTI.A,LTI.B,weight.Q,weight.R,[]);
    K = -K;
    [H,h]=costgen(predmod,weight,dim);
    x(:,1)=LTI.x0;

    fprintf('Starting regulation MPC simulation R:%2.1f\n',R_list(i))
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
        x(:,k+1)=LTI.A*x_0_MPC + LTI.B*u_rec(:,k);
        
    end
    R_results{i,1} = u_rec;
    R_results{i,2} = x(:,1:end-1);
    R_results{i,3} = t;
    clc
    fprintf('Simulation finished\n')


end

figure('name','Comparison of value of R for MPC regulation')
for i=1:length(R_list)
    subplot(2,1,1)
    stairs(R_results{i,3},R_results{i,2}(2,:))
    hold on
    subplot(2,1,2)
    stairs(R_results{i,3},R_results{i,1}(1,:))
    hold on
end

subplot(2,1,1)
xlabel('Time (s)')
ylabel('Position (m)')
grid on
legend(cellstr(num2str(R_list', 'R = %-1.0e')))
subplot(2,1,2)
xlabel('Time (s)')
ylabel('Collective Command (% of maximum)')
grid on
legend(cellstr(num2str(R_list', 'R = %-1.0e')))