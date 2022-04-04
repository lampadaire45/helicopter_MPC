%% Init
clear all
close all
clc

set(0, 'DefaultLineLineWidth', 1.0);

%% Get dynamics
% u = [theta_0 theta_c]
% x = [x z u w q theta lambda_i]
init_hover_dynamics

%% MPC reference
% x = [x   z   u   w   q   theta lambda_i]
x_0 = [0   0   0   0   0   0     0     ]';

% u = [theta_0 theta_c]
u_0 = [0       0]';

%Definition of the LTI system
LTI.A=sysd.A; 
LTI.B=sysd.B;
LTI.C=[eye(2) zeros(2,5)];
LTI.Cd=0;
LTI.x0=x_0;

%Definition of system dimension
dim.nx=size(LTI.A,1);     %state  dimension
dim.nu=size(LTI.B,2);     %input  dimension
dim.ny=size(LTI.C,1);     %output dimension
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

% Generation of prediction model 
predmod=predmodgen(LTI,dim);            
[H,h]=costgen_tracking(predmod,weight,dim);

% Constraints
%   u = [theta_0 theta_c] % percentage of max value
u_lim = [0.25    1]';
%   x = [x    z    u   w   q            theta          lambda_i]
x_lim = [1000 1000 5   5   deg2rad(5)   deg2rad(15)    1000]';
x_lim_vec = repmat(x_lim,[dim.N+1,1]);

% Time vectors
t_span = [0 10];
t = [t_span(1):sysd.Ts:t_span(2)]';

T=length(t);                  %simulation horizon
x=zeros(dim.nx,T+1);
u_rec=zeros(dim.nu,T);
x(:,1)=LTI.x0;

% Reference
r = [zeros(1,(T-1)/2), 0*ones(1,(T-1)/2+1);
     zeros(1,(T-1)/2), 1*ones(1,(T-1)/2+1)];
LTI.last_yref = [rand(size(r,1),1)];

% Receding horizon implementation
for k=1:T
    fprintf('Time is %2.2f \n',t(k))

    x_0=x(:,k);

    % Calculate OTS
    LTI.yref = r(:,k);
    if any(LTI.last_yref ~= LTI.yref)
        [xr,ur]=optimalss(LTI,dim,0);
    end
    LTI.last_yref = r(:,k);

    % Solve problem with CVX
    cvx_begin quiet
        variable uN(dim.nu*dim.N)
        minimize(0.5*uN'*H*uN+(h*[x_0; xr; ur])'*uN)
        
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
    
end

% Plots
plot_single(t,x(:,1:end-1)',u_rec','MPC constrained')

%% Trajectory following (effect of change of Nref)

Nref_list = [1,2,5,10,15,25];
x_Nref = {};
u_Nref = {};
for m=1:length(Nref_list)
dim.Nref = Nref_list(m);   %reference advance knowing

% x = [x   z   u   w   q   theta lambda_i]
x_0 = [0   0   0   0   0   0     0     ]';

% u = [theta_0 theta_c]
u_0 = [0       0]';

%Definition of the LTI system
LTI.A=sysd.A; 
LTI.B=sysd.B;
LTI.C=[eye(2) zeros(2,5)];
LTI.Cd=0;
LTI.x0=x_0;

%Definition of system dimension
dim.nx=size(LTI.A,1);     %state  dimension
dim.nu=size(LTI.B,2);     %input  dimension
dim.ny=size(LTI.C,1);     %output dimension
dim.N=25;       %horizon

%Definition of quadratic cost function
%           x = [    x    z    u     w    q    theta   lambda_i]
weight.Q = 1E2*diag([1E0  1E0  1E0   1E0  1E0  1E0     1E-4]);
%           u = [theta_0 theta_c]
weight.R = diag([1       1]);
weight.beta = 1;

[K,weight.P,e] = dlqr(LTI.A,LTI.B,weight.Q,weight.R,[]);
K = -K;

% Generation of prediction model 
predmod=predmodgen(LTI,dim);            
[H,h]=costgen_tracking(predmod,weight,dim);

% Constraints
%   u = [theta_0 theta_c] % percentage of max value
u_lim = [0.25    1]';
%   x = [x    z    u   w   q            theta          lambda_i]
x_lim = [1000 1000 5   5   deg2rad(5)   deg2rad(15)    1000]';
x_lim_vec = repmat(x_lim,[dim.N+1,1]);

% Time vectors
t_span = [0 10];
t = [t_span(1):sysd.Ts:t_span(2)]';

T=length(t);                  %simulation horizon
x=zeros(dim.nx,T+1);
u_rec=zeros(dim.nu,T);
x(:,1)=LTI.x0;

% Reference
r = [zeros(1,5/sysd.Ts), 0*ones(1,T-(5/sysd.Ts));
     zeros(1,5/sysd.Ts), 0.5*ones(1,T-(5/sysd.Ts))];
LTI.last_yref = [rand(size(r,1),1)];


    % Receding horizon implementation
    for k=1:T
        fprintf('Time is %2.2f \n',t(k))
    
        x_0=x(:,k);
        
        % Select next Nref yref values. Fill if needed
        if k>T-dim.Nref+1
            LTI.yref = r(:,k:T);
            for i=1:k-(T-dim.Nref+1)
                LTI.yref(:,end+1) = LTI.yref(:,end);
            end
        else
            LTI.yref = r(:,k:k+dim.Nref-1);
        end
    
        % Calculate OTS
    %     LTI.yref = r(:,k);
    %     if any(LTI.last_yref ~= LTI.yref)
    %         [xr,ur]=optimalss(LTI,dim,0);
    %     end
    %     LTI.last_yref = r(:,k);
        xr = [];
        ur = [];
        for i=1:dim.Nref
            xr((i-1)*dim.nx+1:i*dim.nx,1) = [LTI.yref(:,i); zeros(5,1)];
            ur((i-1)*dim.nu+1:i*dim.nu,1) = zeros(dim.nu,1);
        end
    
        % Solve problem with CVX
        cvx_begin quiet
            variable uN(dim.nu*dim.N)
            minimize(0.5*uN'*H*uN+(h*[x_0; xr/dim.Nref; ur])'*uN)
            
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
        
    end

x_Nref{m} = x;
u_Nref{m} = u_rec;

end

% Plots
figure
for m=1:length(Nref_list)

    subplot(2,1,1)
    plot(t,-x_Nref{m}(2,1:end-1))
    hold on
    subplot(2,1,2)
    plot(t,u_Nref{m}(1,:))
    hold on
    

end
subplot(2,1,1)
plot(t,-r(2,:),'--')
xlabel('Time (s)')
ylabel('Position (m)')
grid on
legend([cellstr(num2str(Nref_list', 'Nref = %-d'));'Reference'])
axis([0 10 -inf inf])
subplot(2,1,2)
xlabel('Time (s)')
ylabel('Command (% of maximum)')
grid on
legend(cellstr(num2str(Nref_list', 'Nref = %-d')))
axis([0 10 -inf inf])

%% Trajectory following Nref fixed

% x = [x   z   u   w   q   theta lambda_i]
x_0 = [0   0   0   0   0   0     0     ]';

% u = [theta_0 theta_c]
u_0 = [0       0]';

%Definition of the LTI system
LTI.A=sysd.A; 
LTI.B=sysd.B;
LTI.C=[eye(2) zeros(2,5)];
LTI.Cd=0;
LTI.x0=x_0;

%Definition of system dimension
dim.nx=size(LTI.A,1);     %state  dimension
dim.nu=size(LTI.B,2);     %input  dimension
dim.ny=size(LTI.C,1);     %output dimension
dim.N=25;       %horizon
dim.Nref = 5;   %reference advance knowing

%Definition of quadratic cost function
%           x = [    x    z    u     w    q    theta   lambda_i]
weight.Q = 1E2*diag([1E0  1E0  1E0   1E0  1E0  1E0     1E-4]);
%           u = [theta_0 theta_c]
weight.R = diag([1       1]);
weight.beta = 1;

[K,weight.P,e] = dlqr(LTI.A,LTI.B,weight.Q,weight.R,[]);
K = -K;

% Generation of prediction model 
predmod=predmodgen(LTI,dim);            
[H,h]=costgen_tracking(predmod,weight,dim);

% Constraints
%   u = [theta_0 theta_c] % percentage of max value
u_lim = [0.25    1]';
%   x = [x    z    u   w   q            theta          lambda_i]
x_lim = [1000 1000 5   5   deg2rad(5)   deg2rad(15)    1000]';
x_lim_vec = repmat(x_lim,[dim.N+1,1]);

% Time vectors
t_span = [0 15];
t = [t_span(1):sysd.Ts:t_span(2)]';

T=length(t);                  %simulation horizon
x=zeros(dim.nx,T+1);
u_rec=zeros(dim.nu,T);
x(:,1)=LTI.x0;

% Reference
r = [zeros(1,5/sysd.Ts), 5*ones(1,T-(5/sysd.Ts));
     zeros(1,1/sysd.Ts), 10*ones(1,T-(1/sysd.Ts))];
LTI.last_yref = [rand(size(r,1),1)];


    % Receding horizon implementation
    for k=1:T
        fprintf('Time is %2.2f \n',t(k))
    
        x_0=x(:,k);
        
        % Select next Nref yref values. Fill if needed
        if k>T-dim.Nref+1
            LTI.yref = r(:,k:T);
            for i=1:k-(T-dim.Nref+1)
                LTI.yref(:,end+1) = LTI.yref(:,end);
            end
        else
            LTI.yref = r(:,k:k+dim.Nref-1);
        end
    
        % Calculate OTS
    %     LTI.yref = r(:,k);
    %     if any(LTI.last_yref ~= LTI.yref)
    %         [xr,ur]=optimalss(LTI,dim,0);
    %     end
    %     LTI.last_yref = r(:,k);
        xr = [];
        ur = [];
        for i=1:dim.Nref
            xr((i-1)*dim.nx+1:i*dim.nx,1) = [LTI.yref(:,i); zeros(5,1)];
            ur((i-1)*dim.nu+1:i*dim.nu,1) = zeros(dim.nu,1);
        end
    
        % Solve problem with CVX
        cvx_begin quiet
            variable uN(dim.nu*dim.N)
            minimize(0.5*uN'*H*uN+(h*[x_0; xr/dim.Nref; ur])'*uN)
            
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
        
    end

% Plots
plot_position_angle_control(t,x(:,1:end-1)',u_rec','Trajectory following',r')
