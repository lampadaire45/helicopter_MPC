function [admis_x,k,av] = plotXf(Xf,entries,range,states,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if ~isempty(varargin)
    plot_data = varargin{1};
else
    plot_data = true;
end

admis_x = [];
x = zeros(7,1);

for i=1:length(range)
    for j=1:length(range)
        x(entries(1)) = range(i);
        x(entries(2)) = range(j);
        
        if all(Xf.A*x<=Xf.b)
            admis_x(:,end+1) = x;
        end
    end
end

if length(admis_x) == length(range)^2
    warning('Range is most probably to small to contain all true values')
end

[k,av] = convhull(admis_x(entries(1),:),admis_x(entries(2),:));

if plot_data
    figure_name = sprintf('Varying %s and %s while keeping other states at 0',states{entries(1)},states{entries(2)});
    figure('Name',figure_name)
    %plot(admis_x(entries(1),:),admis_x(entries(2),:),'*g')
    hold on
    fill(admis_x(entries(1),k),admis_x(entries(2),k),'g','FaceAlpha',0.3)
    plot(admis_x(entries(1),k),admis_x(entries(2),k),'r','linewidth',2)
    %axis([min(range) max(range) min(range) max(range)])
    axis padded
    xlabel(states{entries(1)})
    ylabel(states{entries(2)})
    terminal_set_text = sprintf('Terminal Set (Area %2.1f)',av);
    legend(terminal_set_text,'Limit of Terminal Set')
end
end