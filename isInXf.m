function [out] = isInXf(x)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
load('Xf.mat');
out = all(Xf.A*x<=Xf.b);
end