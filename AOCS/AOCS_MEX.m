%% AOCS system sizing
clc; close all; clear;

%% Mars express example

% Reaction wheels
h = 12; % Nms
Trw = 0.075; % Nm

% Thrusters
Isp = 292; % s
L = 0.69; % force arm, m
F = 10; % N

% Gravity Gradient data
% Shape: 1.5 m x 1.8 m x 1.1 m
V = 1.5*1.8*1.1; % m^3
Imax = 1088.31; % kg*m^2
Imin = 373.08; % kg*m^2
th_max = 68*pi/180; % max deviation of z axis from radial direction to the planet, rad
mu_mars = 4.28e4; % km^3/s^2

% Orbit
a = 9392; % km
e = 0.61;
b = a * sqrt(1 - e^2); % mean radius, km
R = sqrt(a * b); % A_ellipse = A_circle -> equivalent circle radius, km

% Solar Radiation Pressure data
cp_cg = 0.3; % distance between satellite center of solar pressure and of gravity, m 
Fs = 590; % Sun heat flux at Mars, W/m^2
Asp = 12; % m^2
q = 0.3;

%% Disturbance torques

% GG
T_GG = (3 * mu_mars / (2 * R^3)) * (Imax - Imin) * sin (2 * th_max); % Nm

% SRP
I = 0; % worst case scenario incidence angle
c = 300000000; % speed of light, m/s
T_SRP = (Fs / c) * Asp * (1 + q) * cos(I) * (cp_cg); % Nm

% For MetOp-A we have to add Drag and Magnetic Field

% margin for preliminary estimation / statistical = 100%
T_dist = 2 * (T_GG + T_SRP);

%% Actuators: reaction wheels

% Disturbances rejection
T_op = 2 * pi * sqrt(a^3 / mu_mars); % orbital period, s
h_dist_period = T_dist * T_op; % Nms
desat = h / h_dist_period; % number of orbit after which desaturation is needed

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

% Prpoellant mass
g0 = 9.81;
m_prop_1rw = t_desat_imposed * F_des / (Isp * g0); % to desaturate 1 RW, kg

deltat_des = 72 * 3600; % desaturate each 3 orbits, s
lifetime = 2 * 365 * 24 * 60 * 60; % s
N_RW = 4;
m_prop_des = N_RW * m_prop_1rw * lifetime / deltat_des; % total propellant mass to desaturate, kg


% For MetOp-A we have to calculate also magnetorquers







