function [x_dot_vect] = dynamics(t,states,input,param)
%Dynamics Calculate states derivates of helicopter non-linear system

% u = [theta_0 theta_c]
% x = [x z u w q theta lambda_i]

%% Parameters
rho_alt = param.rho_alt;
g = param.g;
aircraft = param.aircraft;
weight = aircraft.weight;

%% Inputs
if isa(input,'function_handle')
    input = input(t);
end

theta_0 = input(1)*aircraft.collective_max;
theta_c = input(2)*aircraft.cyclic_max;


%% States
x = states(1);
z = states(2);
u = states(3);
w = states(4);
q = states(5);
theta = states(6);
lambda_i = states(7);

tau=0.1;		%time constant in dynamiCs inflow!!!

%% Other parameters
V = sqrt(u^2+w^2);
V_bar = V/(aircraft.omega*aircraft.rotor_r);
if V==0 
    V=0.001;
end
phi = atan2(w,u);
alpha_c = -phi+theta_c;
mu = V_bar*cos(alpha_c);
lambda_c = V_bar*sin(alpha_c);

gamma = atan2(w,u);

a_1 = ((-16*q)/(aircraft.lock*aircraft.omega)+(8/3)*mu*theta_0-2*mu*(lambda_c+lambda_i))/(1-0.5*mu^2);
alpha_d = alpha_c-a_1;

C_T_elem = (aircraft.cla*aircraft.sigma/4)*((2/3)*theta_0*(1+(3/2)*mu^2)-(lambda_c+lambda_i));
C_T_gl = 2*lambda_i*sqrt((V_bar*cos(alpha_d))^2+(V_bar*sin(alpha_d)+lambda_i)^2);

T = C_T_elem*rho_alt*(aircraft.omega*aircraft.rotor_r)^2*pi*aircraft.rotor_r^2;
D_fus = 0.5*rho_alt*aircraft.ACd_fuse*V^2;

%% Forces
Fx = -weight*sin(theta)-D_fus*cos(theta)+T*sin(theta_c-a_1);
Fz = weight*cos(theta)-D_fus*sin(theta)-T*cos(theta_c-a_1);
My = -T*aircraft.h*sin(theta_c-a_1);

%% Derivatives
x_dot_vect = zeros(length(states),1);

x_dot = u*cos(theta)+w*sin(theta);
z_dot = -u*sin(theta)+w*cos(theta);
u_dot = -q*w+Fx/(weight/g);
w_dot = q*u+Fz/(weight/g);
q_dot = My/aircraft.Iyy;
theta_dot = q;
lambda_i_dot = (C_T_elem-C_T_gl)/tau;

x_dot_vect(1) = x_dot;
x_dot_vect(2) = z_dot;
x_dot_vect(3) = u_dot;
x_dot_vect(4) = w_dot;
x_dot_vect(5) = q_dot;
x_dot_vect(6) = theta_dot;
x_dot_vect(7) = lambda_i_dot;

end