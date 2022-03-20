

%% Physical
R = 0.77; % rotor radius (m)
%R_1 = 0.415/2; % Radius (rotor center at beginning pallet) (m)
y_m = 0; % Distance center of mass to main rotor y-axis (m)
l_m = 0; % Distance center of mass to main rotor x-axis (m)
h_m = 0; % Distance center of mass to main rotor z-axis (m)
h_t = 0; % Distance center of mass to tail rotor z-axis (m)
l_t = 0.6; % Distance center of mass to tail rotor x-axis (m)
I = [0.24 0.58 1.6]; % Inertia matrix (kg/m^2)
%omega = 150; % Rotor angular velocity (rad/s)

%rho = 1.225; % Air density (kg/m^3)
%a = 6; % Polar slope (1/rad)
%e = 0; % Hinge (m)
%c = 0.07; % Chord (m)
%B = 2; % Number of blades
m = 7; % Helicopter mass (kg)
g = 9.81; % Gravity (m/s^2)
%M_b = ((m*R^2)/2)*(1-e/R)^2*g; % Static moment (N*m)
%I_b = ((m*R^2)/3)*(1-e/R)^3*g; % Inertia moment (kg*m^2)
%e_MR = 2.25/12; % Main rotor hinge offset (m) CHECK THIS
%gamma = rho*a*c*R^4/I_b;

A_QMR = 0.000018; %Coefficient of ratio
B_QMR = 0.01; % Initial drag of main rotor
