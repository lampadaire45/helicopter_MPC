function [A,B] = lin_sys(x,u,param)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% u = [theta_0 theta_c]
% x = [x z u w q theta lambda_i]
t =  0;
h = 0.01;

A = zeros(length(x),length(x));
B = zeros(length(x),length(u));

for i=1:length(x)
    xp = x; xp(i) = xp(i)+h;
    xm = x; xm(i) = xm(i)-h;

    A(:,i) = (dynamics(t,xp,u,param)-dynamics(t,xm,u,param))/(2*h);

end

for j=1:length(u)
    up = u; up(j) = up(j)+h;
    um = u; um(j) = um(j)-h;

    B(:,j) = (dynamics(t,x,up,param)-dynamics(t,x,um,param))/(2*h);

end

end