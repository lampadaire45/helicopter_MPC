%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Offset free tracking controller
%
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

fprintf('Initialized linear and non-linear system dynamics\n')

%% LTI system definition
% x = [x   z   u   w   q   theta lambda_i]'
x_0 = [0   0   0   0   0   0     0     ]';

% u = [theta_0 theta_c]'
u_0 = [0       0]';

% d = [delta_u delta_w]'
% y = [x y q w]'

%Definition of the LTI system
LTI.A=sysd.A; 
LTI.B=sysd.B;
LTI.C=[eye(2) zeros(2,5);
       0 0 0 0 1 0 0];
LTI.x0=x_0;
LTI.Bd=[eye(2);
        zeros(5,2)];
LTI.Cdopt=[1 0;
         0 1;
         0 0];
LTI.Cd = zeros(3,2);
LTI.d = [0;0];

%Definition of system dimension
dim.nx=size(LTI.A,1);     %state  dimension
dim.nu=size(LTI.B,2);     %input  dimension
dim.ny=size(LTI.C,1);     %output dimension
dim.nd=size(LTI.Cd,2);
dim.N=25;      %horizon
dim.Nref = 1;

%Definition of quadratic cost function
%           x = [    x    z    u     w    q    theta   lambda_i]
weight.Q = 1E2*diag([1E0  1E0  1E0   1E0  1E0  1E0     1E-4]);
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

fprintf('Defining MPC normal system\n')

%% Check if the problem is well posed

if rank(obsv(LTI.A,LTI.C)) == dim.nx
    fprintf('Original system is observable\n')
end

if rank([eye(dim.nx)-LTI.A -LTI.Bd; LTI.C LTI.Cd]) == dim.nx+dim.nd
    fprintf('Augmented system is observable\n')
end

%% Extended system computation

LTIe.A=[LTI.A LTI.Bd; zeros(dim.nd,dim.nx) eye(dim.nd)];
LTIe.B=[LTI.B; zeros(dim.nd,dim.nu)];
LTIe.C=[LTI.C LTI.Cd];
LTIe.x0=[LTI.x0; LTI.d];

%Definition of system dimension
dime.nx=size(LTIe.A,1);     %state  dimension
dime.nu=size(LTIe.B,2);     %input  dimension
dime.ny=size(LTIe.C,1);     %output dimension
dime.N=dim.N;            %horizon
dime.Nref=dim.Nref;      %reference horizon

%Definition of quadratic cost function
weighte.Q=blkdiag(weight.Q,zeros(dim.nd));            %weight on output
weighte.R=weight.R;                                   %weight on input
weighte.P=blkdiag(weight.P,zeros(dim.nd));            %terminal cost
weighte.beta = weight.beta;

x_lime = [x_lim;1E3;1E3];
x_lime_vec = repmat(x_lime,[dim.N+1,1]);

fprintf('Defining MPC augmented system\n')

%% Offset-free MPC from output (Disturbance recognition)
LTIe.x0 = [0 0 0 0 0 0 0 0 0]';

predmode=predmodgen(LTIe,dime);  
[He,he]=costgen_tracking(predmode,weighte,dime);

L=place(LTIe.A',LTIe.C',linspace(0.35,0.50,9))';      %observer gain 

% Kalman observer gains
% Q_kf = diag([1E0 1E0 1E0]);
% R_kf = diag([1 1 1 1 1 1 1 1 1]);
% 
% [~,Obs_eigvals,Obs_gain] = dare(LTIe.A',LTIe.C',R_kf,Q_kf);
% Obs_gain = Obs_gain';
% L = Obs_gain;

% Time vectors
t_span = [0 10];
t = [t_span(1):sysd.Ts:t_span(2)]';

% Receding horizon implementation
T=length(t);
xe=zeros(dime.nx,T+1);
x=zeros(dim.nx,T+1);
y=zeros(dime.ny,T+1);
u_rec=zeros(dime.nu,T);
xehat=zeros(dime.nx,T+1);

xe(:,1)=LTIe.x0;
x(:,1) = LTI.x0;
xehat(:,1)=[0 0 0 0 0 0 0 0 0]';
y(:,1)=LTIe.C*LTIe.x0;


% Reference
r = [0*ones(1,5/sysd.Ts) 0*ones(1,1+5/sysd.Ts);  %x_ref
     0*ones(1,5/sysd.Ts) 0*ones(1,1+5/sysd.Ts);  %z_ref
     0*ones(1,5/sysd.Ts) 0*ones(1,1+5/sysd.Ts)]; %q_ref
LTI.last_yref = rand(size(r,1),1);

% Disturbance
d = [0*ones(1,1/sysd.Ts)  0.1*ones(1,1+9/sysd.Ts);  %u_d
     0*ones(1,7/sysd.Ts)  -0.01*ones(1,1+3/sysd.Ts)];  %w_d

fprintf('Starting disturbance recognition simulation\n')
for k=1:T
    fprintf('Time is %2.2f \n',t(k))

    xe_0=xe(:,k);%+[zeros(7,1); d(:,k)];
    dhat=xehat(end-dim.nd+1:end,k);
    
    if t(k)>=2
        fprintf('Time is 2 sec\n');
    end

    LTI.yref = r(:,k);
    %Compute optimal ss (online, at every iteration)
    [xr,ur]=optimalss(LTI,dim,dhat); 
    xre=[xr;dhat];
    
    % Solve problem with CVX
    cvx_begin quiet
        variable uN(dim.nu*dim.N)
        minimize(0.5*uN'*He*uN+(he*[xehat(:,k); xre; ur])'*uN)
        
        % input constraints
        uN <=  repmat(u_lim,[dim.N,1]);
        uN >= -repmat(u_lim,[dim.N,1]);
        % state constraints
        predmode.S*uN <= -predmode.T*xehat(:,k) + x_lime_vec;
        predmode.S*uN >= -predmode.T*xehat(:,k) - x_lime_vec; 

    cvx_end   

    % Select the first input only
    u_rec(:,k) =uN(1:1*dim.nu);
    
    x(:,k+1) = LTI.A*x(1:dim.nx,k) + LTI.B*u_rec(:,k) + LTI.Bd*d(:,k);
    y(:,k) = LTI.C*x(:,k) + LTI.Cd*d(:,k);

    % Update extended-state estimation
    xehat(:,k+1)=LTIe.A*xehat(:,k)+LTIe.B*u_rec(:,k)+L*(y(:,k)-LTIe.C*xehat(:,k));

end
fprintf('Finished disturbance recognition simulation\n')

% Plots
plot_position_angle_control(t,x(:,1:end-1)',u_rec','Offset Free following')

figure('Name','Tracking error')
e=y(:,1:end-1)-r;
plot(t,e),
xlabel('Time (s)'), ylabel('Tracking error'), grid on;
legend('x','z','q');

figure('name','Disturbance recognition')
stairs(t,xehat(end-1:end,1:end-1)')
hold on
stairs(t,d','--')
xlabel('Time (s)')
ylabel('Disturbance (m)')
legend('x_d estimated','z_d estimated','x_d real','z_d real')
grid on

%% Offset-free MPC from output (Following reference with disturbance)
LTIe.x0 = [0 0 0 0 0 0 0 0 0]';

predmode=predmodgen(LTIe,dime);  
[He,he]=costgen_tracking(predmode,weighte,dime);

L=place(LTIe.A',LTIe.C',linspace(0.35,0.50,9))';      %observer gain 

% Kalman observer gains
% Q_kf = diag([1E0 1E0 1E0]);
% R_kf = diag([1 1 1 1 1 1 1 1 1]);
% 
% [~,Obs_eigvals,Obs_gain] = dare(LTIe.A',LTIe.C',R_kf,Q_kf);
% Obs_gain = Obs_gain';
% L = Obs_gain;

% Time vectors
t_span = [0 20];
t = [t_span(1):sysd.Ts:t_span(2)]';

% Receding horizon implementation
T=length(t);
xe=zeros(dime.nx,T+1);
x=zeros(dim.nx,T+1);
y=zeros(dime.ny,T+1);
u_rec=zeros(dime.nu,T);
xehat=zeros(dime.nx,T+1);

xe(:,1)=LTIe.x0;
x(:,1) = LTI.x0;
xehat(:,1)=[0 0 0 0 0 0 0 0 0]';
y(:,1)=LTIe.C*LTIe.x0;


% Reference
r = [0*ones(1,5/sysd.Ts) 1*ones(1,5/sysd.Ts) 0*ones(1,1+10/sysd.Ts);  %x_ref
     0*ones(1,5/sysd.Ts) 0*ones(1,5/sysd.Ts) 0*ones(1,1+10/sysd.Ts);  %z_ref
     0*ones(1,5/sysd.Ts) 0*ones(1,5/sysd.Ts) 0*ones(1,1+10/sysd.Ts)]; %q_ref
LTI.last_yref = rand(size(r,1),1);

% Disturbance
d = [0*ones(1,1/sysd.Ts)  0.1*ones(1,4/sysd.Ts) 0*ones(1,1+15/sysd.Ts);  %u_d
     0*ones(1,7/sysd.Ts)  0*-0.01*ones(1,3/sysd.Ts) 0*0.1*ones(1,1+10/sysd.Ts)];  %w_d

fprintf('Starting tracking with disturbance simulation\n')
for k=1:T
    fprintf('Time is %2.2f \n',t(k))

    xe_0=xe(:,k);%+[zeros(7,1); d(:,k)];
    dhat=xehat(end-dim.nd+1:end,k);
    
    if t(k)>=2
        fprintf('Time is 2 sec\n');
    end

    LTI.yref = r(:,k);
    %Compute optimal ss (online, at every iteration)
    [xr,ur]=optimalss(LTI,dim,dhat); 
    xre=[xr;dhat];
    
    % Solve problem with CVX
    cvx_begin quiet
        variable uN(dim.nu*dim.N)
        minimize(0.5*uN'*He*uN+(he*[xehat(:,k); xre; ur])'*uN)
        
        % input constraints
        uN <=  repmat(u_lim,[dim.N,1]);
        uN >= -repmat(u_lim,[dim.N,1]);
        % state constraints
        predmode.S*uN <= -predmode.T*xehat(:,k) + x_lime_vec;
        predmode.S*uN >= -predmode.T*xehat(:,k) - x_lime_vec; 

    cvx_end   

    % Select the first input only
    u_rec(:,k) =uN(1:1*dim.nu);
    
    x(:,k+1) = LTI.A*x(1:dim.nx,k) + LTI.B*u_rec(:,k) + LTI.Bd*d(:,k);
    y(:,k) = LTI.C*x(:,k) + LTI.Cd*d(:,k);

    % Update extended-state estimation
    xehat(:,k+1)=LTIe.A*xehat(:,k)+LTIe.B*u_rec(:,k)+L*(y(:,k)-LTIe.C*xehat(:,k));

end
fprintf('Finished tracking with disturbance simulation\n')

% Plots
plot_position_control(t,x(:,1:end-1)',u_rec','Offset Free following',r')

figure('Name','Tracking error')
e=y(:,1:end-1)-r;
plot(t,e),
xlabel('Time (s)'), ylabel('Tracking error'), grid on;
legend('x','z','q');

figure('name','Disturbance recognition')
stairs(t,xehat(end-1:end,1:end-1)')
hold on
stairs(t,d','--')
xlabel('Time (s)')
ylabel('Disturbance (m)')
legend('x_d estimated','z_d estimated','x_d real','z_d real')
grid on