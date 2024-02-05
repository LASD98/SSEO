%% MetOp-A EPS system sizing

clc; close all; clearvars;

%% Variables
SC_mat = 0; % 0 = Si, 1 = Triple Junction Solar Cells InGaP/GaAs/Ge
BAT_mat = 0; % 0 = Li-ion, 1 = Ni-Cd

%% Variable Initialization
% max power requirements of system
P_day = 1810; % W 
P_eclipse = 1810; % W

% voltage
V_sys = 28; % V

% time info
T_day = 67*60; % s
T_eclipse = 34*60; % s
T_life = 5; % years

% line efficiencies
x_d = 0.8; % in daylight
x_e = 0.6; % in eclipse

% irradiance
I0 = 1361; % irradiance of Earth, W/m^2
theta = 30; % inclination angle, worst case 30 deg

switch SC_mat
    case 0
        % performance of ISS silicon solar cells
        E_bol = 0.145; % beginning of life efficiency
        dpy = 0.028; % degradation/year
        I_d = 0.7; % inherent degredation
        % t = 0.05; % thickness, cm
        % rho = 75; % density, kg/m^3
        avw = 320; % average weight, g/m^2
        A_cell = 0.0064; % cell area, m^2
        V_sa_cell = 0.628; % cell voltage, V
    case 1
        % performance of CESI CTJ30 Triple Junction Solar Cells InGaP/GaAs/Ge
        E_bol = 0.295; % beginning of life efficiency
        dpy = 0.0375; % degradation/year
        I_d = 0.7; % inherent degredation
        avw = 890; % average weight, g/m^2
        A_cell = 0.003015; % cell area, m^2
        V_sa_cell = 2.61; % cell voltage, V
end

% battery sizing info (source: slides tables)
N_bat = 5; % number of batteries
C_bat = 40; % battery capacity, Ah

switch BAT_mat
    case 0
        % Li-ion
        DOD = 0.25; % depth of discharge (5 years in LEO)
        eta_bat = 0.4; % transmission efficiency between batteries and load
        E_m = 140; % Wh/kg specific energy
        E_v = 250; % Wh/dm^3 energy density
        V_bat_cell = 3.6; % cell voltage, V
        mu = 0.8; % package efficiency
    case 1
        % Ni-Cd
        DOD = 0.2; % depth of discharge (5 years in LEO)
        eta_bat = 0.4; % transmission efficiency between batteries and load
        E_m = 40; % Wh/kg specific energy
        E_v = 90; % Wh/dm^3 energy density
        V_bat_cell = 1.35; % cell voltage, V
        mu = 0.8; % package efficiency
end

%% Power request to solar array
% in max power demand condition 
P_sa = ((P_eclipse * T_eclipse) / (x_e * T_day)) + (P_day / x_d)

%% Specific Power
P0 = E_bol * I0; % specific power output
P_bol = P0 * I_d * cosd(theta); % SA specific power

L_life = (1 - dpy)^T_life;
P_eol = L_life * P_bol;

%% Physical Measurements of SA
A_sa = P_sa / P_eol;
% m_sa = A_sa * rho_si * t

N_cells = ceil(A_sa / A_cell);
N_sa_series = ceil(V_sys / V_sa_cell);
N_real = ceil(N_cells / N_sa_series) * N_sa_series;

A_sa_real = N_real * A_cell
m_cells = A_sa_real * avw / 1000

%% Battery sizing
% T_r and P_r taken for eclipse
P_r = P_eclipse; % W
T_r = T_eclipse / 3600; % h
C = (T_r * P_r) / (DOD * N_bat * eta_bat);

m_bat = C / E_m
V_bat = C / E_v

N_bat_series = ceil(V_sys / V_bat_cell);
V_real = N_bat_series * V_bat_cell;
C_string = mu * C_bat * V_real;

N_parallel = ceil(C/C_string);
C_real = N_parallel * C_string;

