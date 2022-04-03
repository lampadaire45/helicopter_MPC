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
weight.Q = 1E2*diag([1E0  1E0  1E0   1E0  1E0  1E0     1E0]);
%           u = [theta_0 theta_c]
weight.R = diag([1       1]);

[K,weight.P,e] = dlqr(sysd.A,sysd.B,weight.Q,weight.R,[]);
K = -K;

eig(LTI.A+LTI.B*K);

u_lim = 1*ones(size(LTI.B,2),1);
x_lim = 10*ones(size(LTI.A,1),1);

[Xf, Z] = findXf(LTI.A, LTI.B, K, -x_lim, x_lim, -u_lim, u_lim);

%% Plot Xf
states = {'x (m)','z (m)','u (m/s)', 'w (m/s)','q (rad/s)', '\theta (rad)','\lambda_i'};
range = [-5:0.1:5];

entries = [1 2];
plotXf(Xf,entries,range,states)

entries = [3 4];
plotXf(Xf,entries,range,states)

entries = [5 6];
plotXf(Xf,entries,range,states)

entries = [1 3];
plotXf(Xf,entries,range,states)

entries = [2 4];
plotXf(Xf,entries,range,states)

%% Xn
