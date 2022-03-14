function [A,B] = lin_sys(x,u)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
t =  0;
h = 0.01;

A = zeros(length(x),length(x));
B = zeros(length(x),length(u));

for i=1:length(x)
    xp = x; xp(i) = xp(i)+h;
    xm = x; xm(i) = xm(i)-h;

    A(:,i) = (sys(t,xp,u)-sys(t,xm,u))/(2*h);

end

for j=1:length(u)
    up = u; up(j) = up(j)+h;
    um = u; um(j) = um(j)-h;

    B(:,j) = (sys(t,x,up)-sys(t,x,um))/(2*h);

end

end