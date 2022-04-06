function [u_sat] = sat_input(u)
%sat_input Saturates input
%   u: vector (mx1)

%% Saturation limits
T_MR_max = 0.25; %deviation from hover value
T_TR_max = 0.25; %deviation from hover value
beta_1s_max = 1; %deviation from hover value
beta_1c_max = 1; %deviation from hover value

%% Saturate
u_sat = zeros(size(u));

for i=1:size(u,2)
    u_sat(1,i) = min(max(u(1,i), -T_MR_max), T_MR_max);
    u_sat(2,i) = min(max(u(2,i), -T_TR_max), T_TR_max);
    u_sat(3,i) = min(max(u(3,i), -beta_1s_max), beta_1s_max);
    u_sat(4,i) = min(max(u(4,i), -beta_1c_max), beta_1c_max);
end


end