function [X_dot] = sys(t,states,input)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% x = [x y z u  v  w  phi  theta  psi  p  q  r]
% u = [T_MR  T_TR  beta_1s  beta_1c]

%% Param
param

%% States and inputs
x = states(1);
y = states(2);
z = states(3);
u = states(4);
v = states(5);
w = states(6);
phi = states(7);
theta = states(8);
psi = states(9);
p = states(10);
q = states(11);
r = states(12);

T_MR = input(1);
T_TR = input(2);
beta_1s = input(3);
beta_1c = input(4);

%% System dynamics
x_dot = u;
y_dot = v;
z_dot = w;

u_dot = (-T_MR*sin(beta_1c)-sin(theta)*m*g)/m+v*r-w*q;
v_dot = (T_MR*sin(beta_1s)+sin(phi)*cos(theta)*m*g)/m+u*r+w*p;
w_dot = (-T_MR*cos(beta_1s)*cos(beta_1c)+cos(phi)*cos(theta)*m*g)/m+u*q-v*p;

phi_dot = p+sin(phi)*tan(theta)*q+cos(phi)*tan(theta)*r;
theta_dot = cos(phi)*q-sin(psi)*r;
psi_dot = (sin(phi)/cos(theta))*q+(cos(phi)/cos(theta))*r;

p_dot = ((I(2)-I(3))*q*r+T_MR*cos(beta_1s)*cos(beta_1c)*y_m+T_TR*h_t+T_MR*sin(beta_1s)*h_m-(A_QMR*T_MR^2+B_QMR)*sin(beta_1c))/I(1);
q_dot = ((I(3)-I(1))*p*r+T_MR*cos(beta_1s)*cos(beta_1c)*l_m+T_MR*sin(beta_1c)*h_m-(A_QMR*T_MR^2+B_QMR)*sin(beta_1s))/I(2);
r_dot = ((I(1)-I(2))*p*q+T_MR*sin(beta_1s)*l_m-T_TR*l_t-T_MR*sin(beta_1c)*y_m-(A_QMR*T_MR^2+B_QMR)*cos(beta_1c)*cos(beta_1s))/I(3);

%% Output
X_dot = zeros(12,1);

X_dot(1) = x_dot;
X_dot(2) = y_dot;
X_dot(3) = z_dot;
X_dot(4) = u_dot;
X_dot(5) = v_dot;
X_dot(6) = w_dot;
X_dot(7) = phi_dot;
X_dot(8) = theta_dot;
X_dot(9) = psi_dot;
X_dot(10) = p_dot;
X_dot(11) = q_dot;
X_dot(12) = r_dot;

end