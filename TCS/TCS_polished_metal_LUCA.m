clear all
close all
clc

%% MetOp-A
e = 0.0194701; % orbit eccentricity, https://www.heavens-above.com/orbit.aspx?satid=29499
l1 = 6.3; % m
l2 = 2.5; % m
l3 = 2.5; % m
R_pl = 6371; % km
R_med = R_pl + 824; % km
p = R_med*(1+e*cos(45*pi/180)); % km
Rp = p/(1+e); % pericenter, km
Ra = p/(1-e); % apocenter, km

T_min = 273+15; % min s/c temperature, 15 degreed margin, K
Q_int_min = 1500; % internal power generated, cold case, W
T_max = 314-15; % max s/c temperature, 15 degrees margin, K
Q_int_max = 2000; % internal power generated, hot case, W
K_E = 1; % Earth diffusion factor

% equivalent sphere
A_tot = 2*(l1*l2 + l1*l3 + l2*l3);
r_sphere = sqrt(A_tot/(4*pi));

% Solar flux
q0 = 1367.5; % W/m^2, solar flux at 1 AU (Earth)
q_sun = q0; % distance sc-sun = distance earth-sun

% Albedo
alb = 0.35; % Earth albedo [0.31-0.39], https://ttu-ir.tdl.org/bitstream/handle/2346/72957/ICES_2017_142.pdf
theta = 0; % sunsynchronous orbit (?)
q_alb_max = q_sun*alb*cos(theta)*(R_pl/Rp)^2;

% Infrared 
sigma = 5.67e-8; % W/m^2K^4
eps_E = 0.9; % https://www.jpl.nasa.gov/images/pia18833-nasa-spacecraft-maps-earths-global-emissivity
T_pl = 13.9 + 0.86 + 273.15; % K, https://www.climate.gov/news-features/understanding-climate/climate-change-global-temperature
q_IR_max = sigma*eps_E*T_pl^4*(R_pl/Rp)^2;
q_IR_min = sigma*eps_E*T_pl^4*(R_pl/Ra)^2;

%% HOT CASE
% the choice of eps and alpha are based on considering that the mostly used
% material in Met-Op is aluminum and consulting
% 'MetOp-SG SCA Instrument Thermal Design and Verification'
% that brings to choose Kapton (good thermal properties) 
% so consulting the chart of absorbivity/emissivity materials on the slides
eps = 0.6; % emissivity
alpha =  0.4; % absorbivity
h_min = Rp-R_pl; % Min altitude, km
F_max = 0.5*(1-sqrt((h_min/R_pl)^2+2*h_min/R_pl)/(1+h_min/R_pl)); % view factor

% Q sun
A_cross = pi*r_sphere^2;
Q_sun = A_cross*alpha*q_sun; %W

% Q albedo
Q_alb_max = A_tot*F_max*alpha*K_E*q_alb_max;

% Q infrared
Q_IR_max = A_tot*F_max*q_IR_max;

% Tsc from heat balance Qemitted=Qint+Qsun+Qalb+Qirr
T_sc_hot = ((Q_int_max+Q_sun+Q_alb_max+Q_IR_max)/(sigma*eps*A_tot))^(1/4) % s/c temperature, K
if T_sc_hot > T_max
    disp('T_sc_hot > T_max, temperature outside the range, a cooler must be introduced')
end

% Radiators sizing
eps_rad = 0.88; % rediators emissivity
Q_tot_max = Q_int_max+Q_sun+Q_alb_max+Q_IR_max;
A_rad_min = (Q_tot_max-sigma*eps*A_tot*T_max^4)/(sigma*(eps_rad-eps)*T_max^4) % radiators surface, m^2
A_e = A_tot-A_rad_min; % area not covered by radiators

%% COLD CASE
h_max = Ra-R_pl; % Max altitude, km
F_min = 0.5*(1-sqrt((h_max/R_pl)^2+2*h_max/R_pl)/(1+h_max/R_pl)); % Min view factor

Q_IR_min = A_tot*F_min*q_IR_min;

T_sc_cold = ((Q_int_min+Q_IR_min)/(sigma*(eps*A_tot)))^(1/4) % s/c temperature, K
if T_sc_cold < T_min
    disp('T_sc < T_min, a heater must be introduced')
end

% heaters sizing
Q_heaters = sigma*(eps*A_tot)*T_min^4 - Q_IR_min - Q_int_min 
% too high