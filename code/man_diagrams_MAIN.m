% MANOEUVRING DIAGRAM CODE

clc;
clear all;
close all;


%% DATA DEFINITION

MTOW = 26691; % [kg]
S_w = 70.2; % [m^2]
C_Lcr = 0.33;
C_LmaxTO = 1.8;
C_Lmax = C_LmaxTO - 0.7343;
C_Ncr =  1.1*C_Lcr; 
C_Nmax =  1.1*C_Lmax; 
C_NmaxTO =  1.1*C_LmaxTO; 

R = 8.31432/28.964420e-3; % [N*m*kg^-1*K^-1]
gamma = 1.4;
T_SL=288.15; % [k]
T_cr = 216.65; % [k]
rho_SL=1.225; % [kg/m^3]
rho_cr=0.3108; % [kg/m^3]
Mach_cruise = 0.8;
V_cruise = 236.06; % [m/s]
V_cruise_SL = Mach_cruise*sqrt(gamma*R*T_SL); % [m/s]
V_stall = 58.2; % [m/s]


%% MANOEUVRING DIAGRAM PLOTS

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% SL conditions:
i=2;
man_diagr_density(rho_SL,MTOW,Mach_cruise,S_w,C_Nmax,C_NmaxTO,V_cruise_SL,T_SL,R,gamma,i);

% Cruise conditions:
i=3;
man_diagr_density(rho_SL,MTOW,Mach_cruise,S_w,C_Nmax,C_NmaxTO,V_cruise,T_cr,R,gamma,i);