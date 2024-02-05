%% Propulsion system sizing
clc; close all; clearvars;

%%
flag_pres_gas = 0;  %0 for using He; 1 for N2
flag_material = 1;  %0 for Ti; 1 for Al
flag_shape = 1;     %0 for spherical; 1 for cylindrical

m_dry = 3769;           %[Kg]

m_fin = 1.2 * m_dry;

Dv = 166 + 0.3 * 166;  % Dv + 20 % of margin  [m/s]

% Specific impulse 
I = 230;            %[s]

% Density of Hydrazine
rho = 1.01 * 1e3;         %[Kg/m^3]

g = 9.81;   %[m/s^2];


MR = exp(Dv/(I * g));  %Mass ratio  [-]

m_in = MR * m_fin;

m_prop = m_in - m_fin;

m_prop = m_prop + 0.055 * m_prop;       % 5.5% of margin


%Compute the volume of Hydrazine
V_prop = m_prop/rho;          %[Kg]

V_prop = V_prop + 0.1 * V_prop;     %Adding 10%


if flag_pres_gas == 0
    R = 2077.3;         %[J/Kg K]
    gamma = 1.67;
end

if flag_pres_gas == 1
    R = 296.8;
    gamma = 1.40;
end

%Assuming Blow Down (4-6)
B = 5;          %B = P gas iniz/ P gas fin
Dp_feed = 0.05;         %MPa
P_chamb = 1;  %MPa    
DP_inj = 0.3 * P_chamb;
P_tank_f = P_chamb + DP_inj + Dp_feed;
P_gas_f = P_tank_f;

V_gas_in = V_prop/(B-1);  %Isothermal transformation
V_gas_fin = B * V_gas_in;
P_gas_in = B * P_gas_f;

T_tank = 293;       %[K]

P_tank_in = B * P_tank_f;
V_tank = V_gas_in + V_prop;
V_tank = V_tank / 4;   %No of Tanks in reality
V_tank = V_tank + 0.01 * V_tank;         % 1% of margin for bladder

m_press = (P_tank_in * 1e6 * V_gas_in)/(R * T_tank);

m_press = m_press + 0.2 * m_press;          % 20% of margin

switch flag_material 

    case 0
        rho_tank = 2780;    %[Kg/m^3]
        sigma_tank = 950;   %[MPa]

        case 1
         rho_tank = 2810;    %[Kg/m^3]
         sigma_tank = 503;   %[MPa]
end

        switch flag_shape

            case 0
                r_tank  = (3 * V_tank/(4 * pi))^(1/3);
                t_tank = (P_tank_in * r_tank)/(2 * sigma_tank);
                m_tank = rho_tank * (4/3) * pi * ((r_tank + t_tank)^3 - r_tank^3);

            case 1
               h = linspace(0.1,2,100);      %Length of Cyl tank from 10 cm to 2 m 

               % h = 1;
               r_tank = sqrt(V_tank./(pi .* h));
               t_tank =  (P_tank_in .* r_tank)/( sigma_tank);
               m_tank = rho_tank .* pi .* h .* ((r_tank + t_tank).^2 - r_tank .^2); 

        end
    
        %Thrusters MONARC 22-6
        m_thr = 0.72 * 16;

   m_PS = m_tank * 4 + m_press + m_thr;
   m_PS = m_PS + 0.1 * m_PS;

%% PLOTS
if flag_shape == 1
    
figure(1)
subplot(1,3,1)
plot(h,r_tank);
title('Radius variation');

subplot(1,3,2)
plot(h,m_tank);
title('Mass variation')

subplot(1,3,3)
plot(h,t_tank);
title('Thickness variation');
end





