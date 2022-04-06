function [u] = fu2u(input_u,t,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if isa(input_u,'function_handle')
        for i=1:length(t)
        u(i,:) = input_u(t(i),x(i,:)')';    
        end
    else
        u= input_u;
    end
end