% Code for computing and graphing the Design Point (1st and 2nd iteration)
% The solutions provided are: required wing surface and needed thrust at take-off

% 2020, Aircraft Design

% Authors: 
% Cristian Asensio García
% Juan Garrido Moreno
% Yi Qiang Ji Zhang
% Alexis Leon Delgado
% Alba Molina Cuadrado
% David Morante Torra
% Teresa Peña Mercadé
% Ferran Rubio Vallhonrat
% Iván Sermanoukian Molina
% Santiago Villarroya Calavia

% PREAMBLE 
clc
close all
clear all

% DATA INPUT
% S = 51.1; % Wing surface (1st iteration) [m^2]
S = 70.2; % Wing surface (2nd iteration) [m^2]

% Wing_span=20.9; % Wing span  (1st iteration) [m]
Wing_span=23; % Wing span (2nd iteration) [m]

% MTOW=28643; % MTOW (1st iteration) [kg]
MTOW=26720; % MTOW (2nd iteration) [kg] computed by Multiple Height Methodology with Runge-Kuta

xlimsup=4500; % Upper limit of the x-axis
W_S_ratio=linspace(0,xlimsup,1000); % Wing loading vector

% 1. TAKE-OFF LIMITATION
T_W_ratio_to=8.37E-5*W_S_ratio;
figure
plot(W_S_ratio,T_W_ratio_to,'b','DisplayName','Take-off'); 
hold on

% 2. SECOND SEGMENT LIMITATION
Ne=2; % Number of engines
T_ratio=1;
W2_W_to_ratio=0.98;
DeltaCD_0=0.015;
% A = 8.57; % Aspect ratio (1st iteration)
A = 7.5; % Aspect ratio (2nd iteration)
phi = 0.85;
k=1/(pi*A*phi);
[CD_0_sec,k_sec] = Polar_parameters_function(MTOW,S,Wing_span);
CGR=0.024; % Climb gradient rate
CL_max_to=1.9;
CD=CD_0_sec+DeltaCD_0+k*CL_max_to^2;
T_W_ratio_second=Ne/(Ne-1)*T_ratio*W2_W_to_ratio*(CD/CL_max_to+CGR);
yline(T_W_ratio_second,'m','DisplayName','Second segment');
hold on

% 3. CRUISE LIMITATION
lambda = 4.3;
h = 12; 
Tcr_Tto = (0.0013*lambda-0.0397)*h-0.0248*lambda+0.7125;
Tto_Tcr = Tcr_Tto^-1;
rho = 0.3108;
V_cr = 236.06;
Wcr_Wto = 0.9605; % Weight in cruise vs take-off ratio
W_cr=MTOW*Wcr_Wto;
[Cd0_cr,k_cr] = Polar_parameters_function(MTOW,S,Wing_span);
T_W_ratio_cruise =Tto_Tcr/2*rho*V_cr^2./W_S_ratio.*(Cd0_cr+ (W_S_ratio*Wcr_Wto).^2/((0.5*rho*V_cr^2)^2*pi*A*phi));
plot(W_S_ratio,T_W_ratio_cruise,'g','DisplayName','Cruise');
hold on

% 4. LANDING LIMITATION
% s_l=1750; % Landing distance (1st iteration)
s_l=1670; % Landing distance (2nd iteration)
rho=1.225;
% CL_max_l=2.9; % (1st iteration)
CL_max_l=2.83; % (2nd iteration)
v_A=sqrt(s_l*3.2808/0.6/0.3)*0.51444; %[m/s]
v_sl=v_A/1.3;
Wl_S_ratio=v_sl^2*rho*CL_max_l/2;
Wto_Wl_ratio=0.456;
Wto_S_ratio=Wl_S_ratio*Wto_Wl_ratio;
xline(Wto_S_ratio,'r','DisplayName','Landing');
legend('location','northwest','interpreter','latex');

% 5. DESIGN POINT SELECTION
% Wto_S_ratio_design=4000; % (1st iteration)
Wto_S_ratio_design=3735; % (2nd iteration)
S_w=MTOW*9.81/Wto_S_ratio_design % Required wing surface [m^2]
% Tto_Wto_design=0.436; % (1st iteration)
Tto_Wto_design=0.4675; % (2nd iteration)
Tto=Tto_Wto_design*MTOW*9.81 % Required total thrust [N]

hold on
plot(Wto_S_ratio_design,Tto_Wto_design,'o','MarkerFaceColor','y','DisplayName','Design point');

% Axis format
set(gca,'TickLabelInterpreter','latex','fontsize',10)
xlabel('$W_{to}/S_W\;\left[\mathrm{N}/\mathrm{m^2}\right]$','interpreter','latex','FontSize',12)
ylabel('$T_{to}/W_{to}$','interpreter','latex','FontSize',12)
xlim([0 xlimsup])
ylim([0 1])

% Grid format
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.2;