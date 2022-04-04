%% Physical Parameters
g = 9.81; %m/s^2

aircraft.rotor_r = 11/2; % Rotor radius (m)
Nr = 383; %RPM
aircraft.n_blades = 4; % Number of blades
aircraft.c_blade = 0.325; %Chord for each blade (m)

aircraft.omega = Nr/60*2*pi; %rad/s
aircraft.lock = 6; % arbitrarly, for hingeless system
aircraft.sigma = aircraft.n_blades*aircraft.c_blade/(pi*aircraft.rotor_r);

aircraft.cla=5.7; % CL_alpha (set arbitrarly for NACA 0012)
aircraft.ACd_fuse = 3.75; % Cd_fuselage*frontal area, determined graphically from graph lecture notes P12 W4_L08_P2

aircraft.h = 1.26; %m Distance from rotor to CoG (m) Taken as tailboom height
aircraft.Ixx = 360*g; %kg m^2 Taken from Pamflet HH77
aircraft.Iyy = 400*g; %kg m^2 Taken from Pamflet HH77
aircraft.Izz = 1500*g; %kg m^2 Taken from Pamflet HH77

aircraft.weight = 3585*g; %N
rho_alt = 1.225; %kg/m^3

aircraft.cyclic_max = deg2rad(15);
aircraft.collective_max = deg2rad(20);

%% Set as param
heli_param.rho_alt = rho_alt;
heli_param.g = g;
heli_param.aircraft = aircraft;

%% Calculate initial inflow
% u = [theta_0 theta_c]
% x = [x z u w q theta lambda_i]

u_0 = [0;0];
x_0 = zeros(7,1);


V = x_0(3);
% Hover induced flow
vi_h = sqrt(aircraft.weight/(2*rho_alt*pi*aircraft.rotor_r^2)); %induced velocity (m/s)
% Non-dimensional Flight speed
V_bar = V./vi_h;
% Non-dimensional Induced velocity
vi_bar = sqrt((-V_bar.^2+sqrt(V_bar.^4+4))/2);
% Induced speed
vi = vi_bar*vi_h;
% Non-dimensional Induced velocity
lambda_i = vi./(aircraft.omega*aircraft.rotor_r);
x_0(7) = lambda_i;


%% Trim
% Calculating initial cyclic and collective position
[u_0,fval] = trim(x_0,u_0,heli_param);

%% Linearize
[A,B] = lin_sys(x_0,u_0,heli_param);
sysc = ss(A,B,eye(size(A)),[]);

%% Discretize
h = 0.1;
sysd = c2d(sysc,h);

%% Clear variables
x_0_eq = x_0;
u_0_eq = u_0;
clearvars -except sysc sysd u_0_eq x_0_eq heli_param
