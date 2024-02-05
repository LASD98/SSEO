%% TTMTC system sizing METOP-A
clc; close all; clearvars;


%% Frequency domains

%S-band (TT&C uplink and downlink)
Sul = 2053.4; %uplink, MHz
Sdl = 2230; %downlink, MHz
R_Sul = 2; %data rate uplink, kbit/s
R_Sdl = 4; %data rate downlink, kbit/s

%X-band (global data) 
Xdl = [7750,7900]; %downlink, MHz
R_Xdl = 70 * 10^3; %data rate downlink, kbit/s

%VHF-band (LRPT)
Vdl = 137.1; %downlink, MHz
R_Vdl = 72; %data rate downlink, kbit/s

%L-band (AHRPT)
Ldl = 1701.3; %downlink, MHz
R_Ldl = 3500; %data rate downlink, kbit/s

% Worst Case Scenario Downlink Totals
R_BPSK = R_Sdl % kbit/s
R_QPSK = R_Xdl + R_Vdl + R_Ldl %kbit/s

%% Encoding and modulation

%Encoding : Reed-Solomon (for all data types)
Eb_N0_min = 5.5; % minimum Eb/N0 for R-S with BER=10^-5, dB
a_enc = 1.14; 

%Modulation
a_BPSK = 1; %For TT&C
a_QPSK = 2; %For global data, LRPT and AHRPT

% bit error rate
BER_TM_dl = 1e-5; % BER for TM and data downlink
BER_TC_ul = 1e-7; % BER for TC uplink

%% Real data rate (** not sure about this **) 
R_real = (R_BPSK * a_enc / a_BPSK) + (R_QOSK * a_enc / a_QPSK) %kbit/s



%% %%% Equations from MEX %%% %%

% transmitted power
mu_amp = 0.56; % efficiency of the amplifier (TWTA at 90 W, depends on input power)
Ptx = mu_amp * XP;
Ptx_dbW = 10*log10(Ptx); % power in dbW


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

% error per bit to noise density ratio
Eb_N0 = P_rx-N0-10*log10(R)
Eb_N0_min = Eb_N0_min + 3 % Eb/N0 min + 3 dB margin

disp('Eb_N0 must be higher than Eb_N0_min')


%% carrier modulation index reduction
B_mod = 78*pi/180; % between 45 and 90 deg, usually 78
% best efficiency with B_mod = 90 deg,
% normal-mode downlink bit rates -> B_mod = 60 deg,
% Lower B_mod for very low bit rates (e.g. emergency mode links).
P_mod_loss = 20*log10(cos(B_mod)); % dB
P_carrier = P_rx+P_mod_loss; % carrier power, dB

% signal to noise ratio SNR
B = 30; % receiver bandwidth (depends on the modulation, usually around 30 Hz)
SNR_carrier = P_carrier - N0 - 10*log10(B); % dB
SNR_min = 10; % minimum SNR depending on the ground station constraints (DSN->10 dB)
SNR_margin = SNR_carrier-SNR_min % dB

disp('SNR_margin must be higher than 3 dB')

