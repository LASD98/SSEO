clc; close all; clearvars;

% total instruments data rate, bps
IDR = (7.5 + 2.1 + 60 + 622 + 40 + 60 + 2.9 + 1500 + 3.9 + 2.4 + 0.166) * 1000;

% storage, bit
ST = 24e9;

% downlink data rate, bit/s
DDR = 70e6;

% orbit duration, s
T = 101.4 * 60;

% time of data accumulation and time od communication needed to downlink, s
rap = DDR/IDR;
t_com = T/rap

% data accumulated and downlinked in one orbit, bit
d_acc = IDR * T
d_dl = DDR * t_com 

% telemtry downlink and telecommand and ranging uplink, bit

TM = 4096 * t_com
TCeR = 2000 * t_com
