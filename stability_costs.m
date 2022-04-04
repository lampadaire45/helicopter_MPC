function [Vf,l] = stability_costs(x,u,weight)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

Vf = zeros(1,length(x)); %Final cost
l = zeros(1,length(x)); % Stage cost

for i=1:length(x)
    Vf(i) = x(:,i)'*weight.P*x(:,i);
    l(i) = x(:,i)'*weight.Q*x(:,i)+u(:,i)'*weight.R*u(:,i);
end

end