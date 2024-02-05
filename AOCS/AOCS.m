%% AOCS system sizing
clc; close all; clear;

%% Metop

% Reaction wheels
h = 40; % Nms
Trw = 0.24; % Nm

% Thrusters
Isp = 230; % s
L = 1; % force arm, m
F = 23.5; % N

% Gravity Gradient data
% Shape: 1.5 m x 1.8 m x 1.1 m
V = 17.6*6.7*5.4; % m^3
Imax = 2088.31; % kg*m^2
Imin = 773.08; % kg*m^2
th_max = 68*pi/180; % max deviation of z axis from radial direction to the planet, rad
mu = 398600; % km^3/s^2

% Orbit
R_e = 6371;
T = 101.4 * 60;  %s
a = (T^2 * mu/(4 * pi^2))^(1/3);
e = 0;
b = a * sqrt(1 - e^2); % mean radius, km
R = sqrt(a * b); % A_ellipse = A_circle -> equivalent circle radius, km
lifetime = 4 * 365 * 24 * 60 * 60; % s

% Solar Radiation Pressure data
cp_cg = 17.6/3; % distance between satellite center of solar pressure and of gravity, m 
Fs = 1373; % Sun heat flux at Mars, W/m^2
Asp = 12; % m^2
q = 0.3;

%% Disturbance torques

% GG
T_GG = (3 * mu / (2 * R^3)) * (Imax - Imin) * sin (2 * th_max); % Nm

% SRP
I = 0; % worst case scenario incidence angle
c = 300000000; % speed of light, m/s
T_SRP = (Fs / c) * Asp * (1 + q) * cos(I) * (cp_cg); % Nm

% Magnetic
M = 7.96*1e15;      %Tm^3, half of earth magn moment on the equator
R = (824+6378)*1e3; %m
B = 2*M/(R^3);      %earth magn field
D = 5;              %Am^2, residual dipole on s/c, range[1-20]
T_magn = D*B;       %Nm

%Drag
rho = 1.13*1e-14;  %kg/m^3, from 09_attitudeslides_chart
cd = 2.2;
mu = 398600*1e9;   %m^3/s^2
A_drag = 6.3*2.5;  %assuming dimensions of launch config
v = sqrt(mu/R);    %approx circular
cpa_cg = 2.5/2;
T_drag = 1/2*rho*cd*A_drag*(v^2)*cpa_cg; %Nm, drag is the same for each surface 

% margin for preliminary estimation / statistical = 100%
Tdist = 2 * (T_GG + T_SRP + T_drag + T_magn);

%% Actuators: reaction wheels

% Disturbances rejection
T_op = T; % orbital period, s
h_dist_period = Tdist * T_op; % Nms
desat = h / h_dist_period; % number of orbit after which desaturation is needed
n_orb = lifetime/T; % n_ orb during lifetime
n_desat = n_orb/desat; %n of time where desat occurr during lifetime

% Slew maneuver of th_m = 180° around max axis of inertia, worst case
th_m = 180; 
th_dot_max = 0.5; % max slew rate
t = th_m / th_dot_max; 
T_slew = 4 * th_m*pi/180 * Imax / t^2;
% > Trw -> can't slew with this rate

% max affordable slew rate
t_min = sqrt(4 * th_m*pi/180 * Imax / Trw);
th_dot_max_ach = th_m / t_min; % deg/s

%% Thrusters

% Reaction Wheels desaturation
n = 2; % number
t_burn = h / (n * L * F); % s
t_desat = t_burn; % s

% We prefer higher desaturation time
t_desat_imposed = 5; % s
F_des = h / (n * L * t_desat_imposed); % thrust to desaturate in 5 seconds, N

% Propellant mass
g0 = 9.81;
m_prop_1rw = t_desat_imposed * F_des / (Isp * g0); % to desaturate 1 RW, kg

deltat_des = T_op * desat; % desaturate each desat orbits, s
N_RW = 3;
m_prop_des = N_RW * m_prop_1rw * lifetime / deltat_des; % total propellant mass to desaturate, kg

% Magnetotorques
D_nec = Tdist/B;  %dipole necessary to counteract disturbances







