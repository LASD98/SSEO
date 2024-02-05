%% TCS system sizing
clc; close all; clear;

%% Mars express example

% data
r_Earth = 1; %AU
r_sc = 1.5; %AU
q0 = 1367.5; % solar flux at Earth, W/m^2
theta = 0; % irradiance angle between s/c and planet
Q_int_max = 410; % internal power generated, hot case, W
Q_int_min = 300; % internal power generated, cold case, W
T_max = 308 - 15; % max s/c temperature, 15 degrees margin, K
T_min = 268 + 15; % min s/c temperature, 15 degreed margin, K

% S/C data
l1 = 1.5; % m
l2 = 1.8; % m
l3 = 1.4; % m

% Orbit
a = 9392; % km
e = 0.61;
b = a * sqrt(1 - e^2); % mean radius, km
R = sqrt(a * b); % A_ellipse = A_circle -> equivalent circle radius, km
rp = a*(1-e); % pericenter, km
ra = a*(1+e); % apocenter, km

% sphere model
A_tot = 2*(l1*l2+l2*l3+l3*l1); % total surface, m^2
r_sphere = sqrt(A_tot/(4*pi)); % m

% Thermal heat sources

% Solar flux
q_sun = q0*(r_Earth/r_sc)^2; % W/m^2

% Albedo
al = 0.15; % albedo factor
r_Mars = 3390; % Mars radius, m
% max
q_albedo_max = q_sun*al*cos(theta)*(r_Mars/rp)^2; % W/m^2
% min
q_albedo_min = q_sun*al*cos(theta)*(r_Mars/ra)^2; % W/m^2

% IR
sigma = 5.67*10^-8; % Stefan-Boltzmann constant, W/m^2*k^4
eps_Mars = 0.85; % Mars infrared emissivity
T_Mars = 210; % Mars surface temperature, K
% max
q_IR_max = sigma*eps_Mars*T_Mars^4*(r_Mars/rp)^2;
% min
q_IR_min = sigma*eps_Mars*T_Mars^4*(r_Mars/ra)^2;

%% Hot case
A_cross = pi*r_sphere^2; % cross sectional area, m^2
alpha = 0.2; % s/c solar absorpivity
Q_Sun = A_cross*alpha*q_sun;

K_Mars = 1; % Mars diffusion factor
h_min = rp-r_Mars; % Min altitude, m
F_Mars_sc_max = 0.5*(1-sqrt((h_min/r_Mars)^2+2*h_min/r_Mars)/(1+h_min/r_Mars)); % Max view factor
Q_albedo_max = A_tot*F_Mars_sc_max*alpha*K_Mars*q_albedo_max;

Q_IR_max = A_tot*F_Mars_sc_max*q_IR_max;

eps_sc = 0.15; % spacecraft emissivity
T_sc_hot = ((Q_int_max+Q_Sun+Q_albedo_max+Q_IR_max)/(sigma*eps_sc*A_tot))^(1/4) % s/c temperature, K
if T_sc_hot > T_max
    disp('T_sc_hot > T_max, temperature outside the range, a cooler must be introduced')
end

% Radiators sizing
eps_rad = 0.88; % rediators emissivity
Q_tot_max = Q_int_max+Q_Sun+Q_albedo_max+Q_IR_max;
A_rad_min = (Q_tot_max-sigma*eps_sc*A_tot*T_max^4)/(sigma*(eps_rad-eps_sc)*T_max^4) % radiators surface, m^2
A_e = A_tot-A_rad_min; % area not covered by radiators

%% Cold case
h_max = ra-r_Mars; % Max altitude, m
F_Mars_sc_min = 0.5*(1-sqrt((h_max/r_Mars)^2+2*h_max/r_Mars)/(1+h_max/r_Mars)); % Min view factor
Q_IR_min = A_tot*F_Mars_sc_min*q_IR_min;

T_sc_cold = ((Q_int_min+Q_IR_min)/(sigma*(eps_sc*A_e+eps_rad*A_rad_min)))^(1/4) % s/c temperature, K
if T_sc_cold < T_min
    disp('T_sc < T_min, a heater must be introduced')
end

% heaters sizing
Q_heaters = sigma*(eps_sc*A_e+eps_rad*A_rad_min)*T_min^4 - Q_IR_min - Q_int_min 
% too high, heaters can be introduced just for the more important and low
% temperature sensitive components










