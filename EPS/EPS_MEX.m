%% EPS system sizing

%% MEX Example
clc; close all; clearvars;

%% Variable Initialization
% max power requirements of system
P_day = 500; % W 
P_eclipse = 300; % W

% voltage
V_sys = 28; % V

% time info
T_day = 21480; % s
T_eclipse = 5520; % s
T_life = 2.5; % years

% line efficiencies
x_d = 0.8; % in daylight
x_e = 0.6; % in eclipse

% irradiance
I0 = 588.5; % irradiance of Mars
theta = 30; % inclination angle, worst case 30 deg

% performance of silicon solar cells
E_bol = 0.26; % beginning of life efficiency
dpy = 0.028; % degradation / year
I_d = 0.7; % inherent degredation
t = 0.05; % cm
rho_si = 75; % kg/m^3
A_cell = 0.0035; % m^2
V_sa_cell = 2.6; % V

% battery sizing info
N_bat = 3;
C_cell = 22.5 / 3; % Ah --> divided by 3 to get per battery 
DOD = 0.6; % high end of Lithium battery range
eta_bat = 0.4; % efficiency
E_m = 140; % Wh/kg specific energy
E_v = 250; % Wh/dm^3 energy density
V_bat_cell = 3.6; % V
mu = 0.8; % package efficiency

%% Power request to solar array
% in max power demand condition 
P_sa = ((P_eclipse * T_eclipse) / (x_e * T_day)) + (P_day / x_d);

%% Specific Power
P0 = E_bol * I0; % specific power output
P_bol = P0 * I_d * cosd(theta); % SA specific power

L_life = (1 - dpy)^T_life;
P_eol = L_life * P_bol;

%% Physical Measurements of SA
A_sa = P_sa / P_eol;
m_sa = A_sa * rho_si * t;

N_cells = ceil(A_sa / A_cell);
N_sa_series = ceil(V_sys / V_sa_cell);
N_real = ceil(N_cells / N_sa_series) * N_sa_series;

A_sa_real = N_real * A_cell;

%% Battery sizing
% T_r and P_r taken for eclipse
P_r = P_eclipse; %W
T_r = T_eclipse / 3600; %h
C = (T_r * P_r) / (DOD * N_bat * eta_bat);

m_bat = C / E_m;
V_bat = C / E_v;

N_bat_series = ceil(V_sys / V_bat_cell);
V_real = N_bat_series * V_bat_cell;
C_string = mu * C_cell * V_real;

N_parallel = ceil(C/C_string);
C_real = N_parallel * C_string;

