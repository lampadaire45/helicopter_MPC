function [x_dot] = sys(t,x,input)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% x = [u  v  w  phi  theta  psi  p  q  r]
% u = [T_MR  T_TR  beta_1s  beta_1c]

%% Param
param

%% States and inputs
u = x(1);
v = x(2);
w = x(3);
phi = x(4);
theta = x(5);
psi = x(6);
p = x(7);
q = x(8);
r = x(9);

T_MR = input(1);
T_TR = input(2);
beta_1s = input(3);
beta_1c = input(4);

%% System dynamics
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
x_dot = zeros(9,1);

x_dot(1) = u_dot;
x_dot(2) = v_dot;
x_dot(3) = w_dot;
x_dot(4) = phi_dot;
x_dot(5) = theta_dot;
x_dot(6) = psi_dot;
x_dot(7) = p_dot;
x_dot(8) = q_dot;
x_dot(9) = r_dot;

end