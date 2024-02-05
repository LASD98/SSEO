%% TTMTC system sizing
clc; close all; clearvars;

%% Mars express example

% worst case
TM = 50; % telemtry, kbps
SD = 114; % science data, kbps
R = TM + SD;  

% HGA X-band antenna
Xul = 7.1; % uplink, GHz
Xdl = 8.4; % downlink, GHz
XP = 90; % power, W
Xd = 1.65; % diameter, m
Xm = 25; % mass, kg

% LGA S-band antenna
Sul = 2.1; % uplink, GHz
Sdl = 2.3; % downlink, GHz
SP = 10; % power, W
Sd = 0.4; % diameter, m
Sm = 1.2; % mass, kg

% patch UHF antenna
UHFul = 401; % uplink, MHz
UHFdl = 437; % downlink, MHz

% approximanted telemetry/commands data-rate
TCdl = 2000; % bps
TCul = 500; % bps

% transmitted power
mu_amp = 0.56; % efficiency of the amplifier (TWTA at 90 W, depends on input power)
Ptx = mu_amp * XP;
Ptx_dbW = 10*log10(Ptx); % power in dbW

% bit error rate
BER_TM_dl = 1e-5; % BER for TM and data downlink
BER_TC_ul = 1e-7; % BER for TC uplink

% modulation
a_BPSK = 1;
a_QPSK = 2;

% encoding (Reed-Solomon)
Eb_N0_min = 5.5; % minimum Eb/N0 for R-S with BER=10^-5, dB
a_RS = 1.14;

R_real = R*a_RS/a_QPSK;

% satellite antenna
mu_par = 0.55;
mu_helix = 0.7;
mu_horn = 0.52; % efficiencies

c = 300000000; % speed of light, m/s
lambda = c/(Xdl*1e9); % wavelength, m
G_ant = 10*log10(pi*Xd^2*mu_par/lambda^2); % antenna gain

% ground station (the one we could use, Estrack, is on Ex3 slide 14)
D_ant = 70; % Deep Space Network andtenna diameter, m
G_rx = 68.15; % DSN antenna gain, dB
th_rx = 65.3*lambda/D_ant; % beamwidth, deg

% losses
l_cable = -2; % cable losses, between -1 and -3 dB
r = 2.5e11; % Earth-Mars distance, m
l_space = 20*log10(lambda/(4*pi*r)); % free space losses, dB
eta = 0.01; % pointing accuracy, deg
l_point = -12*(eta/th_rx)^2; % pointing losses, dB
l_atm = -0.05; % atmospheric losses from charts at f=8.4 GHz, dB

% Effective Isotropic Radiated Power
EIRP = Ptx_dbW+G_ant+l_cable; % dB

% received power
P_rx = EIRP+G_rx+l_space+l_atm+l_point; % dB

% System Noise Density
k = 1.38e-23; % Boltzmann constant, Ws/K
Ts = 21; % Temperature, K
N0 = 10*log10(k*Ts); % dB

% error per bit to noise density
Eb_N0 = P_rx-N0-10*log10(R)
Eb_N0_min = Eb_N0_min + 3 % Eb/N0 min + 3 dB margin

% Eb_N0 should be higher than Eb_N0_min

%%









