%% TTMTC system sizing
clc; close all; clearvars;

%% X-band antenna

% information bit-rate
R = 70000000; % 70 Mbps, bps

% X-band antenna
Xdl = 7.8; % downlink, GHz
XP = 90; % power, W
Xd = 0.1; % diameter, m

% transmitted power
mu_amp = 0.56; % efficiency of the amplifier (TWTA at 90 W)
Ptx = mu_amp * XP;
Ptx_dbW = 10*log10(Ptx); % power in dbW

% bit error rate
BER_TM_dl = 1e-5; % BER for TM and data downlink
BER_TC_ul = 1e-7; % BER for TC uplink

% modulation
a_QPSK = 2;

% encoding (Reed-Solomon)
Eb_N0_min = 5.5; % minimum Eb/N0 for R-S with BER=10^-5, dB
a_RS = 1.14;

R = R*a_RS/a_QPSK;

% satellite antenna
mu_par = 0.55; % efficiency of parabolic antenna
c = 300000000; % speed of light, m/s
lambda = c/(Xdl*1e9); % wavelength, m
G_ant = 10*log10(pi*Xd^2*mu_par/lambda^2); % antenna gain

% ground station 
D_ant = 11; % Svalbard ground station X-band antenna diameter, m
G_rx = 10*log10(pi*D_ant^2*mu_par/lambda^2); % ground antenna gain, dB
th_rx = 65.3*lambda/D_ant; % beamwidth, deg

% losses
l_cable = -2; % cable losses, between -1 and -3 dB
r = 2890000; % Slant range, m
l_space = 20*log10(lambda/(4*pi*r)); % free space losses, dB
l_point = -1; % pointing losses, dB
l_atm = -3.85; % atmospheric losses from charts, dB

% Effective Isotropic Radiated Power
EIRP = Ptx_dbW+G_ant+l_cable; % dB

% received power
P_rx = EIRP+G_rx+l_space+l_atm+l_point; % dB

% System Noise Density
k = 1.38e-23; % Boltzmann constant, Ws/K
Ts = 250; % sensor temperature, K
N0 = 10*log10(k*Ts); % dB

% error per bit to noise density ratio
Eb_N0 = P_rx-N0-10*log10(R)
Eb_N0_min = Eb_N0_min + 3 % Eb/N0 min + 3 dB margin

disp('Eb_N0 must be higher than Eb_N0_min')

%% carrier modulation index reduction
B_mod = 60*pi/180; % Modulation index, 60 deg for normal-mode downlink bit rates
P_mod_loss = 20*log10(cos(B_mod)); % dB
P_carrier = P_rx+P_mod_loss; % carrier power, dB

% signal to noise ratio SNR
B = 63000000; % bandwidth, Hz
SNR_carrier = P_carrier - N0 - 10*log10(B); % dB
SNR_min = 10; % minimum SNR depending on the ground station constraints
SNR_margin = SNR_carrier-SNR_min % dB

disp('SNR_margin must be higher than 3 dB')

