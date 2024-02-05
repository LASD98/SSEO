clear all
close all
clc

% MetOpA
e = 0.01; % to distinguish hot and cold cases
l1 = 6.3;
l2 = 2.5;
l3 = 2.5;
R_pl = 6371;
R_med = R_pl + 824;
p = R_med*(1+e*cos(45*3.14/180));
Rp = p/(1+e);
Ra = p/(1-e);

% Data that I don't know ???
T_min = 273.15-20; % K
Q_int_min = 400; % random
T_max = 40+273.15;
Q_int_max = 500; % ???, random value
K_E = 1; % ???,Earth diffusion factor

% equivalent sphere
A_tot = 2*(l1*l2 + l1*l3 + l2*l3);
r_sphere = sqrt(A_tot/(4*pi));

% q sun
q0 = 1367.5; % W/m^2, solar flux at 1AU (Earth)
q_sun = q0; % distance sc-sun = distance earth-sun

% q albedo
alb = 0.35; % earth albedo in range [0.31-0.39]
theta = 0; % sunsynchronous orbit (?)
q_alb_max = q_sun*alb*cos(theta)*(R_pl/Rp)^2;

% q infrared 
sigma = 5.67*1e-8; %W/m^2K^4
eps_E = 0.9; %https://www.jpl.nasa.gov/images/pia18833-nasa-spacecraft-maps-earths-global-emissivity
T_pl = 13.9 + 0.86 + 273.15; %K, https://www.climate.gov/news-features/understanding-climate/climate-change-global-temperature
q_IR_max = sigma*eps_E*T_pl^4*(R_pl/Rp)^2;
q_IR_min = sigma*eps_E*T_pl^4*(R_pl/Ra)^2;

% HOT CASE
% the choice of eps and alpha are based on considering that
% the core of sc is an absorber while external surfaces 
% (made in polished metal mainly) are reflectors
% so consulting the chart of absorbivity/emissivity materials on the slides
eps = 0.1; %emissivity
alpha =  0.2; %absorbivity
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
T_sc_hot = ((Q_int_max+Q_sun+Q_alb_max+Q_IR_max)/(sigma*eps*A_tot))^(1/4); % s/c temperature, K
% T_sc > T_max, temperature outside the range, a cooler must be introduced

% Radiators sizing
eps_rad = 0.85; % radiators emissivity
Q_tot_max = Q_int_max+Q_sun+Q_alb_max+Q_IR_max;
A_rad_min = (Q_tot_max-sigma*eps*A_tot*T_max^4)/(sigma*(eps_rad-eps)*T_max^4);% radiators surface, m^2
A_e = A_tot-A_rad_min; % area not covered by radiators


% COLD CASE
h_max = Ra-R_pl; % Max altitude, km
F_min = 0.5*(1-sqrt((h_max/R_pl)^2+2*h_max/R_pl)/(1+h_max/R_pl)); % Min view factor

Q_IR_min = A_tot*F_min*q_IR_min;

T_sc_cold = ((Q_int_min+Q_IR_min)/(sigma*(eps*A_e+eps_rad*A_rad_min)))^(1/4); % s/c temperature, K
% T_sc > T_min
