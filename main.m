%% Init
clear all
close all
clc

set(0, 'DefaultLineLineWidth', 1.0);
%% Simulation Conditions
% x = [x y z u  v  w  phi  theta  psi  p  q  r]
% u = [T_MR  T_TR  beta_1s  beta_1c]

t_span = [0 10];
x_0 = zeros(12,1);
u_0 = [0.60 -0.1 0 0]';
[u_0,fval] = trim(x_0,u_0);

% Set step function
u_step_col = @(t,x) u_0 +[0.25*heaviside(t-1)-0.25*heaviside(t-2) 0 0 0]';
u_step_cyc = @(t,x) u_0 +[0 0 0.25*heaviside(t-1)-0.25*heaviside(t-2) 0]';

%% Non-linear simulation
f_step_col = @(t,x) sys(t,x,u_step_col);
f_step_cyc = @(t,x) sys(t,x,u_step_cyc);

[t_nl_col,x_nl_col] = ode45(f_step_col,t_span,x_0);
[t_nl_cyc,x_nl_cyc] = ode45(f_step_cyc,t_span,x_0);

%% Plot Non-linear results
pos_axis    = [0 10 -1e1  1e1];
speed_axis  = [0 10 -5e0  5e0];
rate_axis   = [0 10 -5e0  5e0];
angle_axis  = [0 10 -1e1  1e1];
axis = [pos_axis;speed_axis;rate_axis;angle_axis];

plot_single(t_nl_col,x_nl_col,u_step_col,'Non-linear model collective step')
plot_single(t_nl_cyc,x_nl_cyc,u_step_cyc,'Non-linear model cyclic step')
% Highly coupled system

%visualize_helicopter_trajectory_rotating(x_nl_col,[],0.01)

%% Linearize system in hover
% x = [x y z u  v  w  phi  theta  psi  p  q  r]
% u = [T_MR  T_TR  beta_1s  beta_1c]
x_0 = zeros(12,1);
u_0 = [0.60 -0.1 0 0]';
[u_0,fval] = trim(x_0,u_0);

[A,B] = lin_sys(x_0,u_0);
sys_lin = ss(A,B,eye(size(A)),zeros(size(B)));

%% Discretize system
dt = 0.1;

sysd = c2d(sys_lin,dt);

%% Simulate linear
t_lin = [t_span(1):(t_span(2)-t_span(1))/100:t_span(2)]';
[~,t_lin_col,x_lin_col] = lsim(sys_lin,fu2u(u_step_col,t_lin,zeros(length(t_lin),length(x_0)))-u_0'.*ones(length(t_lin),4),t_lin);
[~,t_lin_cyc,x_lin_cyc] = lsim(sys_lin,fu2u(u_step_cyc,t_lin,zeros(length(t_lin),length(x_0)))-u_0'.*ones(length(t_lin),4),t_lin);

%% Plot linear system
plot_single(t_lin_col,x_lin_col,u_step_col,'Linear model collective step')
plot_single(t_lin_cyc,x_lin_cyc,u_step_cyc,'Linear model cyclic step')

%% Compare linear and non-linear system
plot_compare_nl2l(t_lin_col,x_lin_col,u_step_col,t_nl_col,x_nl_col,u_step_col,'Comparison of Linear and Non-Linar System with Collective Step',axis)
plot_compare_nl2l(t_lin_cyc,x_lin_cyc,u_step_cyc,t_nl_cyc,x_nl_cyc,u_step_cyc,'Comparison of Linear and Non-Linar System with Collective Step',axis)

%% Check controlability of system (continuous)

if rank(ctrb(sys_lin.A,sys_lin.B)) == size(sys_lin.A,1)
    fprintf("System is controllable\n")
else
    fprintf("System not controllable\n")
end

%% Check controlability of system (discrete)

if rank(ctrb(sysd.A,sysd.B)) == size(sysd.A,1)
    fprintf("System is controllable\n")
else
    fprintf("System not controllable\n")
end

%% Check stability of system (discrete)
if any(eig(sysd.A)>1)
    fprintf('System is unstable\n')
elseif any(eig(sysd.A)==1)
    fprintf('System marginally stable\n')
else
    fprintf('System is stable\n')
end

%% Plot param
pos_axis    = [0 10 -1.25e0  1.25e0];
speed_axis  = [0 10 -5e0  5e0];
rate_axis   = [0 10 -5e0  5e0];
angle_axis  = [0 10 -5e0  5e0];
axis = [pos_axis;speed_axis;rate_axis;angle_axis];

%% LQR controller
% x = [x y z u  v  w    phi    theta     psi      p  q  r]
x_0 = [1 1 0 0  0  0    0      0         0        0  0  0]';

% u = [T_MR  T_TR  beta_1s  beta_1c]
u_0 = [0     0     0        0]';

t_span = [0 30];
t = [t_span(1):dt:t_span(2)]';
%    x = [x    y    z     u    v    w    phi    theta     psi      p    q    r]
weight.Q = diag([1E2  1E2  1E2   1    1     1   1      1         1        1E1  1E1  1E1]);
%    u = [T_MR  T_TR  beta_1s  beta_1c]
weight.R = diag([1     1     1        1]);

[K,P,e] = dlqr(sysd.A,sysd.B,weight.Q,weight.R,[]);
K = -K;

sysd_LQR_cl = ss(sysd.A+sysd.B*K,sysd.B,sysd.C,sysd.D,dt);

% Check stability
if any(eig(sysd.A+sysd.B*K)>1)
    fprintf('System is unstable\n')
else
    fprintf('System is stable\n')
end

t_LQR = t;
x_LQR = zeros(length(x_0),length(t_LQR));
x_LQR(:,1) = x_0;
u_LQR = zeros(length(u_0),length(t_LQR));

for i=1:length(t)-1
    % Calculate u
    u_LQR(:,i) = K*x_LQR(:,i);

    x_LQR(:,i+1) = sysd.A*x_LQR(:,i)+sysd.B*u_LQR(:,i);

end

plot_single(t_LQR,x_LQR',u_LQR','LQR',axis)

%visualize_helicopter_trajectory_rotating(x_LQR,[],0.5)

%% Constrained LQR

t_LQR_cstrn = t;
x_LQR_cstrn = zeros(length(x_0),length(t_LQR_cstrn));
x_LQR_cstrn(:,1) = x_0;
u_LQR_cstrn = zeros(length(u_0),length(t_LQR_cstrn));

for i=1:length(t)-1
    % Calculate u
    u_LQR_cstrn(:,i) = K*x_LQR_cstrn(:,i);
    
    % Saturate u
    u_LQR_cstrn(:,i) = sat_input(u_LQR_cstrn(:,i));

    x_LQR_cstrn(:,i+1) = sysd.A*x_LQR_cstrn(:,i)+sysd.B*u_LQR_cstrn(:,i);

end

plot_single(t_LQR_cstrn,x_LQR_cstrn',u_LQR_cstrn','LQR_constrained',axis)

%% MPC controller non constrained
% x = [x y z u  v  w    phi    theta     psi      p  q  r]
x_0 = [1 1 0 0  0  0    0      0         0        0  0  0]';

% u = [T_MR  T_TR  beta_1s  beta_1c]
u_0 = [0     0     0        0]';

%Definition of the LTI system
LTI.A=sysd.A; 
LTI.B=sysd.B;
LTI.x0=x_0;

%Definition of system dimension
dim.nx=size(LTI.A,1);     %state dimension
dim.nu=size(LTI.B,2);     %input dimension
dim.N=25;      %horizon

%Definition of quadratic cost function
%    x = [x    y    z     u    v    w    phi    theta     psi      p    q    r]
weight.Q = diag([1E2  1E2  1E2   1    1     1   1      1         1        1E1  1E1  1E1]);
%    u = [T_MR  T_TR  beta_1s  beta_1c]
weight.R = diag([1     1     1        1]);
weight.beta = 3;

% Terminal cost
[K,weight.P,e] = dlqr(sysd.A,sysd.B,weight.Q,weight.R,[]);
K = -K;

% Generation of prediction model 
predmod=predmodgen(LTI,dim);            
[H,h]=costgen(predmod,weight,dim);

% Time vectors
t_span = [0 10];
t = [t_span(1):dt:t_span(2)]';

T=length(t);                  %simulation horizon
x=zeros(dim.nx,T+1);
u_rec=zeros(dim.nu,T);
x(:,1)=LTI.x0;

% Receding horizon implementation
for k=1:T
    fprintf('Time is %2.2f \n',k*dt)

    x_0=x(:,k);
    
%     % Solve the unconstrained optimization problem (with YALMIP)
%     u_uncon = sdpvar(dim.nu*dim.N,1);                        %define optimization variable
%     Constraint=[];                                           %define constraints
%     Objective = 0.5*u_uncon'*H*u_uncon+(h*x_0)'*u_uncon;     %define cost function
%     optimize(Constraint,Objective);                          %solve the problem
%     u_uncon=value(u_uncon);                                  %assign the solution
%     % Select the first input only
%     u_rec(:,k)=u_uncon(1:1*dim.nu);

    % Solve problem with CVX
    cvx_begin quiet
        variable uN(dim.nu*dim.N)
        minimize(0.5*uN'*H*uN+(h*x_0)'*uN)
    cvx_end
    % Select the first input only
    u_rec(:,k)=uN(1:1*dim.nu);

    % Compute the state/output evolution
    x(:,k+1)=LTI.A*x_0 + LTI.B*u_rec(:,k);
    clear u_uncon
    
end

% Plots
plot_single(t,x(:,1:end-1)',u_rec','MPC non-constrained',axis)

%% MPC constrained

% x = [x y z u  v  w    phi    theta     psi      p  q  r]
x_0 = [1 1 0 0  0  0    0      0         0        0  0  0]';

% u = [T_MR  T_TR  beta_1s  beta_1c]
u_0 = [0     0     0        0]';

%Definition of the LTI system
LTI.A=sysd.A; 
LTI.B=sysd.B;
LTI.C=sysd.C;
LTI.x0=x_0;

%Definition of system dimension
dim.nx=size(LTI.A,1);     %state  dimension
dim.nu=size(LTI.B,2);     %input  dimension
dim.ny=size(LTI.C,1);     %output dimension
dim.N=25;      %horizon

%Definition of quadratic cost function
%           x = [x    y    z     u    v    w    phi    theta     psi      p    q    r]
weight.Q = diag([1E0  1E0  1E0   1    1    1   1       1         1        1E0  1E0  1E0]);
%    u = [T_MR  T_TR  beta_1s  beta_1c]
weight.R = diag([1     1     1        1]);
weight.beta = 3;

% Terminal costplot_single(t_LQR,x_LQR',u_LQR','LQR',axis)
[K,weight.P,e] = dlqr(sysd.A,sysd.B,weight.Q,weight.R,[]);
K = -K;

% Generation of prediction model 
predmod=predmodgen(LTI,dim);            
[H,h]=costgen(predmod,weight,dim);

% Constraints
% u =   [T_MR  T_TR  beta_1s  beta_1c]
u_lim = [0.25  0.25  1        1]';
% x =   [x   y   z   u   v  w    phi          theta       psi          p            q          r         ]
x_lim = [1E3 1E3 1E3 5   5  5    deg2rad(10)  deg2rad(10) deg2rad(10)  deg2rad(1)   deg2rad(1) deg2rad(1)]';
x_lim_vec = repmat(x_lim,[dim.N+1,1]);

% Time vectors
t_span = [0 30];
t = [t_span(1):dt:t_span(2)]';

T=length(t);                  %simulation horizon
x=zeros(dim.nx,T+1);
u_rec=zeros(dim.nu,T);
x(:,1)=LTI.x0;

% Receding horizon implementation
for k=1:T
    fprintf('Time is %2.2f \n',k*dt)

    x_0=x(:,k);
    
%     % Solve the unconstrained optimization problem (with YALMIP)
%     u_uncon = sdpvar(dim.nu*dim.N,1);                        %define optimization variable
%     Constraint=[];                                           %define constraints
%     Objective = 0.5*u_uncon'*H*u_uncon+(h*x_0)'*u_uncon;     %define cost function
%     optimize(Constraint,Objective);                          %solve the problem
%     u_uncon=value(u_uncon);                                  %assign the solution
%     % Select the first input only
%     u_rec(:,k)=u_uncon(1:1*dim.nu);

    % Solve problem with CVX
    cvx_begin quiet
        variable uN(dim.nu*dim.N)
        minimize(0.5*uN'*H*uN+(h*x_0)'*uN)
        
        % input constraints
        uN <=  repmat(u_lim,[dim.N,1]);
        uN >= -repmat(u_lim,[dim.N,1]);
        % state constraints
        predmod.S*uN <= -predmod.T*x_0 + x_lim_vec;
        predmod.S*uN >= -predmod.T*x_0 - x_lim_vec; 

    cvx_end
    % Select the first input only
    u_rec(:,k)=uN(1:1*dim.nu);

    % Compute the state/output evolution
    x(:,k+1)=LTI.A*x_0 + LTI.B*u_rec(:,k);
    clear u_uncon
    
end

% Plots
plot_single(t,x(:,1:end-1)',u_rec','MPC non-constrained',axis)

%% Stability
%Definition of the LTI system
LTI.A=sysd.A; 
LTI.B=sysd.B;
LTI.C=sysd.C;
LTI.x0=x_0;

%Definition of quadratic cost function
%           x = [x    y    z     u    v    w    phi    theta     psi      p    q    r]
weight.Q = 1E2*diag([1E0  1E0  1E0   1    1    1   1       1         1        1E0  1E0  1E0]);
%           u = [T_MR  T_TR  beta_1s  beta_1c]
weight.R = diag([1     1     1        1]);

[K,weight.P,e] = dlqr(sysd.A,sysd.B,weight.Q,weight.R,[]);
K = -K;

eig(LTI.A+LTI.B*K)

u_lim = 1*ones(size(LTI.B,2),1);
x_lim = 10*ones(size(LTI.A,1),1);
[A_set,b_set] = max_admissible_set(LTI.A,K,u_lim,x_lim);

[Xf, Z] = findXf(LTI.A, LTI.B, K, -x_lim, x_lim, -u_lim, u_lim);