function [u_trim,fval] = trim(x_0,u_0,param)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%options = optimset('Display','iter','PlotFcns',@optimplotfval);
options = optimset('TolFun',1E-6,'TolX',1E-6);
[u_trim,fval] = fminsearch(@(u)trim_fun(u,x_0,param),u_0,options);


end