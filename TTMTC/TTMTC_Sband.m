%% TTMTC system sizing
clc; close all; clearvars;

%% S-band antenna

Sul = 2053.4583; % uplink, MHz
Sdl = 2230; % downlink, MHz
R_ul = 2000; % data rate uplink, bit/s
R_dl = 4096; % data rate downlink, bit/s
SP = 30; % input power, W
Sd = 0.115; % diameter, m

% Transmitted power 
mu_amp = 0.15; % efficiency of the amplifier (SSA at 30 W)
Ptx = mu_amp*SP; % W
Ptx_dbW = 10*log10(Ptx); % dbW
Ptx_ul = 20000; % ground station transmitted power, W
Ptx_ul_dbW = 10*log10(Ptx_ul); % dbW

% bit error rate
BER_ul = 1e-7; % uplink
BER_dl = 1e-5; % downlink

% modulation
a_QPSK = 2;

% encoding (Reed-Solomon)
Eb_N0_min = 5.5; % minimum Eb/N0 for R-S with BER=10^-5, dB
a_RS = 1.14;

% Real data rates
R_ul = R_ul*a_RS/a_QPSK; % kbit/s
R_dl = R_dl*a_RS/a_QPSK; % kbit/s

% satellite antenna
mu_par = 0.55; % efficiency of parabolic antenna
c = 300000000; % speed of light, m/s
lambda_ul = c/(Sul*1e6); % uplink wavelength, m
lambda_dl = c/(Sdl*1e6); % downlink wavelength, m
G_ant_ul = 10*log10(pi*Sd^2*mu_par/lambda_ul^2); % uplink antenna gain, dB
G_ant_dl = 10*log10(pi*Sd^2*mu_par/lambda_dl^2); % downlink antenna gain, dB

% ground station 
D_ant = 11; % Svalbard ground station S-band antenna diameter, m
G_rx_ul = 10*log10(pi*D_ant^2*mu_par/lambda_ul^2); % ground antenna gain (uplink), dB
G_rx_dl = 10*log10(pi*D_ant^2*mu_par/lambda_dl^2); % ground antenna gain (downlink), dB
th_rx_ul = 65.3*lambda_ul/D_ant; % beamwidth (uplink), deg
th_rx_dl = 65.3*lambda_dl/D_ant; % beamwidth (downlink), deg

% Losses
l_cable = -1; % cable losses, between -1 and -3 dB
r = 2890000; % Slant range, m
l_space_ul = 20*log10(lambda_ul/(4*pi*r)); % free space losses (uplink), dB
l_space_dl = 20*log10(lambda_dl/(4*pi*r)); % free space losses (downlink), dB
eta = 0.01; % pointing accuracy, deg
l_point_ul = -12*(eta/th_rx_ul)^2; %pointing losses (uplink), dB 
l_point_dl = -12*(eta/th_rx_dl)^2; %pointing losses (uplink), dB 
l_atm_ul = -0.4; % atmospheric losses (uplink), dB
l_atm_dl = -0.4; % atmospheric losses (downlink), dB

% Effective Isotropic Radiated Power
EIRP_ul = Ptx_ul_dbW+G_rx_ul+l_cable; % uplink, dB
EIRP_dl = Ptx_dbW+G_ant_dl+l_cable; % uplink, dB

% Received power
P_rx_ul = EIRP_ul+G_ant_ul+l_space_ul+l_atm_ul+l_point_ul; % uplink, dB
P_rx_dl = EIRP_dl+G_rx_dl+l_space_dl+l_atm_dl+l_point_dl; % uplink, dB

% System Noise Density
k = 1.38e-23; % Boltzmann constant, Ws/K
Ts = 250; % sensor temperature, K
N0 = 10*log10(k*Ts); % dB

% Error per bit to noise density ratio
Eb_N0_ul = P_rx_ul-N0-10*log10(R_ul) %uplink
Eb_N0_dl = P_rx_dl-N0-10*log10(R_dl) %downlink
Eb_N0_min = Eb_N0_min + 3 % Eb/N0 min + 3 dB margin

disp('Eb_N0 must be higher than Eb_N0_min')

% carrier modulation index reduction
B_mod_ul = 78*pi/180; % Modulation index (uplink), rad 
B_mod_dl = 60*pi/180; % Modulation index (downlink), rad 
P_mod_loss_ul = 20*log10(cos(B_mod_ul)); % uplink, dB
P_mod_loss_dl = 20*log10(cos(B_mod_dl)); % uplink, dB
P_carrier_ul = P_rx_ul+P_mod_loss_ul; % carrier power (uplink), dB
P_carrier_dl = P_rx_dl+P_mod_loss_dl; % carrier power (downlink), dB

% signal to noise ratio SNR
B_ul = 1500000; % uplink bandwidth (TC&ranging), Hz
B_dl = 2000000; % downlink bandwidth (TM), Hz
SNR_carrier_ul = P_carrier_ul - N0 - 10*log10(B_ul); % uplink, dB
SNR_carrier_dl = P_carrier_dl - N0 - 10*log10(B_dl); % uplink, dB
SNR_min = 10; % minimum SNR depending on the ground station constraints
SNR_margin_ul = SNR_carrier_ul-SNR_min % dB
SNR_margin_dl = SNR_carrier_dl-SNR_min % dB

disp('SNR_margin must be higher than 3 dB')