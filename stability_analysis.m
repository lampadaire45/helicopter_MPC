%% Init
clear all
close all
clc

set(0, 'DefaultLineLineWidth', 1.0);

%% Get dynamics
% u = [theta_0 theta_c]
% x = [x z u w q theta lambda_i]
init_hover_dynamics

%% Stability
%Definition of the LTI system
LTI.A=sysd.A; 
LTI.B=sysd.B;
LTI.C=sysd.C;

%Definition of quadratic cost function
%           x = [    x    z    u     w    q    theta   lambda_i]
weight.Q = 1E2*diag([1E0  1E0  1E0   1E0  1E0  1E0     1E-4]);
%           u = [theta_0 theta_c]
weight.R = diag([1       1]);

[K,weight.P,e] = dlqr(sysd.A,sysd.B,weight.Q,weight.R,[]);
K = -K;

eig(LTI.A+LTI.B*K);

%   u = [theta_0 theta_c] % percentage of max value
u_lim = [0.25    1]';
%   x = [x    z    u   w   q            theta          lambda_i]
x_lim = [1000 1000 5   5   deg2rad(5)   deg2rad(15)    1000]';

[Xf, Z] = findXf(LTI.A, LTI.B, K, -x_lim, x_lim, -u_lim, u_lim);
save('Xf.mat','Xf')
%% Plot Xf
states = {'x (m)','z (m)','u (m/s)', 'w (m/s)','q (rad/s)', '\theta (rad)','\lambda_i'};
range = [-5:0.1:5];

entries = [1 2];
[~,~,~] = plotXf(Xf,entries,range,states);

entries = [3 4];
[~,~,~] = plotXf(Xf,entries,range,states);

entries = [5 6];
[~,~,~] = plotXf(Xf,entries,range,states);

entries = [1 3];
[~,~,~] = plotXf(Xf,entries,range,states);

entries = [2 4];
[~,~,~] = plotXf(Xf,entries,range,states);

%% MPC Constrained

% x = [x   z   u   w   q   theta lambda_i]
x_0 = [1   1 0   0   0   0     0     ]';

% u = [theta_0 theta_c]
u_0 = [0       0]';

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
weight.Q;
weight.R;
weight.beta = 1;

% Terminal costplot_single(t_LQR,x_LQR',u_LQR','LQR',axis)
K;
weight.P;

% Generation of prediction model 
predmod=predmodgen(LTI,dim);            
[H,h]=costgen(predmod,weight,dim);

% Constraints
u_lim;
x_lim;
x_lim_vec = repmat(x_lim,[dim.N+1,1]);

% Time vectors
t_span = [0 10];
t = [t_span(1):sysd.Ts:t_span(2)]';

T=length(t);                  %simulation horizon
x=zeros(dim.nx,T+1);
u_rec=zeros(dim.nu,T);
x(:,1)=LTI.x0;


% Receding horizon implementation
for k=1:T
    fprintf('Time is %2.2f \n',t(k))

    x_0=x(:,k);

%     A_opt = [kron(eye(dim.N),eye(dim.nu));
%             kron(eye(dim.N),-eye(dim.nu));
%             predmod.S;
%             -predmod.S];
%     b_opt = [repmat(u_lim,[dim.N,1]);
%              repmat(u_lim,[dim.N,1]);
%              -predmod.T*x_0 + x_lim_vec;
%               predmod.T*x_0 + x_lim_vec];

    
%     % Solve the unconstrained optimization problem (with YALMIP)
%     u_uncon = sdpvar(dim.nu*dim.N,1);                        %define optimization variable
%     Constraint=[];                                           %define constraints
%     Objective = 0.5*u_uncon'*H*u_uncon+(h*x_0)'*u_uncon;     %define cost function
%     optimize(Constraint,Objective);                          %solve the problem
%     u_uncon=value(u_uncon);                                  %assign the solution
%     % Select the first input only
%     u_rec(:,k)=u_uncon(1:1*dim.nu);
    
    %u_sec = quadprog(H,h*x_0,A_opt,b_opt);

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
plot_single(t,x(:,1:end-1)',u_rec','MPC constrained')

%% Evolution of terminal and stage cost
[Vf,stage_cost] = stability_costs(x(:,1:end-1),u_rec,weight);

figure('name','Evolution of terminal and stage cost')
stairs(t(1:end-1),Vf(2:end)-Vf(1:end-1))
hold on
stairs(t,-stage_cost,'--')
xlabel('Time (s)')
legend('V_f(f(x,u)) - V_f(x)','-l(x,u)')
axis([0 5 -inf inf])

%% Xn calculated for position and keeping all states to zero
% Determined empirically

beta_list = [1 3 10];

x_list = linspace(-1.5,1.5,20);
z_list = linspace(-20,20,20);

xzInXn = zeros(length(x_list), length(z_list), length(beta_list));

for k=1:length(beta_list)
    for i=1:length(x_list)
        for j=1:length(z_list)
            x=zeros(dim.nx,dim.N+1);
            u_rec=zeros(dim.nu,dim.N+1);
            x(:,1)=[x_list(i) z_list(j) zeros(1,5)]';
            for m=1:dim.N
                fprintf('beta: %2.2f x: %2.2f z: %2.2f N is %2.2f \n',beta_list(k), x_list(i), z_list(j),m)
                x_0=x(:,m);
                
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
                u_rec(:,m)=uN(1:1*dim.nu);
            
                % Compute the state/output evolution
                x(:,m+1)=LTI.A*x_0 + LTI.B*u_rec(:,m);

            end
            
            % Check if final value is in Xf
            if all(Xf.A*x(:,end)<=Xf.b)
                xzInXn(i,j,k) = [x_list(i); z_list(j)];
            end


        end
    end
end

save('Xn.mat','xzInXn','x_list','z_list')

%% Plot Xn
load('Xn.mat','xzInXn','x_list','z_list')
x_list;
z_list;
xzInXn;

states = {'x (m)','z (m)','u (m/s)', 'w (m/s)','q (rad/s)', '\theta (rad)','\lambda_i'};
range = [-5:0.1:5];
entries = [1 2];
[Xf_pts,k_Xf,av_Xf] = plotXf(Xf,entries,range,states,false);

figure('Name','Xn for x and z by setting all other states to 0')
[k,av] = convhull(xzInXn(1,:),xzInXn(2,:));
fill(xzInXn(1,k),xzInXn(2,k),'b','FaceAlpha',0.3)
hold on
fill(Xf_pts(1,k_Xf),Xf_pts(2,k_Xf),'g','FaceAlpha',0.3)
axis padded
xlabel(states{entries(1)})
ylabel(states{entries(2)})
terminal_set_text = sprintf('Terminal Set (Area %2.1f)',av_Xf);
region_attraction_text = sprintf('Region of attraction (Area %2.1f)',av);
legend(region_attraction_text,terminal_set_text)
