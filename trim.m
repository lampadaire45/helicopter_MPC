function [u_trim,fval] = trim(x_0,u_0)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%options = optimset('Display','iter','PlotFcns',@optimplotfval);
%[u_trim,fval] = fminsearch(@(u)trim_fun(u,x_0),u_0,options);
[u_trim,fval] = fminsearch(@(u)trim_fun(u,x_0),u_0);


end