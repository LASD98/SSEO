%% MetOp-A OBDH system sizing

clc; close all; clearvars;

C_ADCS = 5*(3*1000 + 16*600 + 2*1500 + 4*500 + 4*800 + 2*1000 + 2000 + 1000 + 15000 + 24000 + 3500 + 13000);
D_ADCS = 5*(3*300 + 16*400 + 2*800 + 4*100 + 4*500 + 2*200 + 200 + 100 + 3500 + 4200 + 2500 + 4000);
TP_ADCS = 5*(3*2.5 + 16*0.6 + 2*6 + 4*1 + 4*4.5 + 2*0.5 + 15 + 12 + 150 + 60 + 4 + 20);

C_PS = 5*(12*800 + 4*800);
D_PS = 5*(12*1500 + 4*1500);
TP_PS = 5*(12*3 + 4*3);

C_EPS = 5*(1200 + 1200);
D_EPS = 5*(500 + 500);
TP_EPS = 5*(0.5 + 0.5);

C_TCS = 5*(400);
D_TCS = 5*(750);
TP_TCS = 5*(1.5);

C_TTMTC = 5*(1000 + 4*1000);
D_TTMTC = 5*(4000 + 4*2500);
TP_TTMTC = 5*(7 + 4*3);

C_OS = 5*(2000 + 700 + 1200 + 3500 + 8000 + 15000 + 4000 + 2000);
D_OS = 5*(700 + 400 + 200 + 2000 + 4000 + 10000 + 1000 + 10000);
TP_OS = 5*(10 + 0.5 + 0.5 + 60 + 60 + 20 + 3 + 5);

C_TOT = C_ADCS + C_PS + C_EPS + C_TCS + C_TTMTC + C_OS
D_TOT = D_ADCS + D_PS + D_EPS + D_TCS + D_TTMTC + D_OS
TP_TOT = TP_ADCS + TP_PS + TP_EPS + TP_TCS + TP_TTMTC + TP_OS

ROM = C_TOT * 16 / 8000
RAM = (C_TOT + D_TOT) * 16 / 8000
TP = TP_TOT / 1000