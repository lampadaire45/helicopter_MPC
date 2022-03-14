function [cost] = trim_fun(u,x_0)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
t = 0;
[x_dot] = sys(t,x_0,u);

cost = x_dot'*x_dot;

end