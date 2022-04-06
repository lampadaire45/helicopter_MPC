function [cost] = trim_fun(u,x_0,param)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% u = [theta_0 theta_c]
% x = [x z u w q theta lambda_i]


t = 0;
[x_dot] = dynamics(t,x_0,u,param);

cost = x_dot(3:end)'*x_dot(3:end); %don't look at derivative of position

end