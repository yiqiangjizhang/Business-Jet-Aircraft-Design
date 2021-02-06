% Code for obtaining the weight design parameters

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

clear
clc
close all
format long

%% DATA INPUT

% Fuselage dimensions
b_f_expected=2.02;
a_f_expected=1.819125;
l_f_expected=21.376375;

% Payload parameters
pax = 10 ; % Number of passengers
W_pax = 77 ; % Mass of a passenger [kg]
W_bagg = 20 ; % Mass of a passenger's baggage [kg]
rho_cargo = 160 ; % Cargo density [kg/m^3]
rho_bagg = 200 ;  % Baggage density [kg/m^3]
kbd = 0.85 ; % Hold occupation efficiency

% Performance data
% Expected total range
total_range=4500*unitsratio('meter','nauticalmile');  % [m]
M_cruise = 0.8; % Cruise Mach number
R_cruise = 0.8*total_range; % Cruise range [m] 
Efficiency = 15 ; % Aerodynamic efficiency during the cruise
h_cruise = 12000 ; % Cruise height [m]
cj_cruise = 2*1.8 * 10^(-5); % Cruise specific consumption [kg/N/s]
cj_loiter = 2*1.1 * 10^(-5); % Loiter specific consumption [kg/N/s]
E_loiter = 30 * 60 ; % Loiter time [s]
E_reserve= 30 * 60 ; % Reserved time to reach alternate aerodrome [s]
MTOW_estimated = 30000 ; % MTOW estimation [kg]
engine_weight=634.225; % Average engine weight [kg]

% Termo-physical properties computation
[T_cruise,rho_cruise]=ISA_atmosphere(h_cruise);
R_g = 2.870528738362446e+02 ; % Constant for air
lambda = 1.4 ; % Heat relation for air
g = 9.81 ; % Gravity acceleration [m/s^2]

% Data regression loading
load('k_b_baggage.mat'); % Constant for the baggage volume
load('alpha.mat'); % Constant for the similarities criterion
load('Delta_W_e_regression.mat'); % Term for the Torenbeek criterion

%% WEIGHT FRACTIONS

% Baggage volume
V_bagg=k_b_baggage*b_f_expected^2*l_f_expected;

% Maximum payload calculation
MPL = pax * (W_pax+W_bagg) + (kbd*V_bagg -((W_bagg*pax)/rho_bagg))*rho_cargo ;  % Computation fo Payload Weight

% Fuel fractions determination for each stage

fraction1 = 0.990 ; % Mass fuel fraction engine start phase
fraction2 = 0.995 ; % Mass fuel fraction taxi phase
fraction3 = 0.995 ; % Mass fuel fraction take off phase
fraction4 = 0.980 ; % Mass fuel fraction Climb phase

% Constant in Breguet's equation for range
ktf = Efficiency * M_cruise*sqrt(lambda*T_cruise*R_g)/(g*cj_cruise);

fraction5_analytical = exp(-R_cruise/ktf); % Mass fuel fraction Cruise phase analytical computation
fraction6_analytical = exp(-E_loiter/(Efficiency/(g*cj_loiter))); % Mass fuel fraction loiter phase analytical computation
fraction_reserve=exp(-E_reserve/(Efficiency/(g*cj_loiter))); % Mass fuel fraction for reaching alternate aerodrome

fraction7 = 0.990 ; % Mass fuel fraction descend phase
fraction8 = 0.992; % Mass fuel fraction shutdown phase

% Mass fuel fraction for the nominal flight
M_ff_analytical = fraction1*fraction2*fraction3*fraction4*fraction5_analytical*fraction6_analytical*fraction7*fraction8*fraction_reserve;  % Total Mass fuel fraction consumed (analytical computation)

%% TORENBEEK CRITERION LOOP FOLLOWING ROSKAM'S METHOD (SINGLE HEIGHT ANALYSIS SHA)

% Delta W_e abscissa calculation
x_fuselage = l_f_expected * 0.5 * (a_f_expected+b_f_expected);
% Delta_e computation
Delta_e = 10^Delta_W_e_regression(2) * x_fuselage^(Delta_W_e_regression(1));

% Initial error definition
error=true;
% Assignation of the estimated MTOW
MTOW_SHA_torenbeek=MTOW_estimated;
% Counter initialization
counter=0;
% Data storage matrix initialization
data_storage_SHA_torenbeek=[];

while error==true % Loop that is executed until convergence is achieved
    
counter=counter+1; % Counter addition

% Total fuel weight computation
W_f_analytical = (1-M_ff_analytical)*MTOW_SHA_torenbeek; 

% Operating empty weight computation (tentative)
OEW_tentative=MTOW_SHA_torenbeek-W_f_analytical-MPL;

% Allowable operating empty weight computation (according to Torenbeek criterion)
OEW_torenbeek = 0.2*MTOW_SHA_torenbeek + 2*engine_weight+ Delta_e+500;

% Relative deviation computation
OEW_torenbeek_deviation=abs(OEW_torenbeek-OEW_tentative)/OEW_torenbeek;

% Data storage
data_storage_SHA_torenbeek(1,counter)=OEW_tentative; % Tentative result
data_storage_SHA_torenbeek(2,counter)=OEW_torenbeek; % Torenbeek result

error=false;

% If the convergence criterion is not fulfilled, the error is set back to
% "true"
if OEW_torenbeek_deviation>0.0005
       error=true;
       
       % Depending on the value, a different adjustment is made
       if OEW_torenbeek>OEW_tentative 
            MTOW_SHA_torenbeek=MTOW_SHA_torenbeek+1;
       else
            MTOW_SHA_torenbeek=MTOW_SHA_torenbeek-1;
       end
       
end % End of the conditional that checks the convergence

end % End of the loop

% Graphic plotting 
fig1=figure(1);
set(fig1,'Renderer', 'painters', 'Position', [400 400 500 350]);
hold on

plot([1:1:size(data_storage_SHA_torenbeek,2)], data_storage_SHA_torenbeek(1,:),'b')

plot([1:1:size(data_storage_SHA_torenbeek,2)], data_storage_SHA_torenbeek(2,:),'r')

% Axis format
set(gca,'TickLabelInterpreter','latex','fontsize',13)
ylabel('$OEW\,\left[\mathrm{kg}\right]$','interpreter','latex','FontSize',15)
xlabel('Number of iterations','interpreter','latex','FontSize',15)
 
% Grid format
grid on
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.2;
grid minor

legend('Tentative','Torenbeek','Location','southeast','interpreter','latex','fontsize',13)

%% SIMILARITIES CRITERION LOOP FOLLOWING ROSKAM'S METHOD (SINGLE HEIGHT ANALYSIS)

% Initial error definition
error=true;
% Assignation of the estimated MTOW
MTOW_SHA_similarities=MTOW_estimated;
% Counter initialization
counter=0;
% Data storage matrix initialization
data_storage_SHA_similarities=[];

while error==true % Loop that is executed until convergence is achieved
    
counter=counter+1; % Counter addition

% Total fuel weight computation
W_f_analytical = (1-M_ff_analytical)*MTOW_SHA_similarities ; 

% Operating empty weight computation (tentative)
OEW_tentative=MTOW_SHA_similarities-W_f_analytical-MPL;

% Allowable operating empty weight (according to similarities criterion)
OEW_similarities = alpha*MTOW_SHA_similarities;

% Relative deviation computation
OEW_similarities_deviation=abs(OEW_similarities-OEW_tentative)/OEW_similarities;

% Data storage
data_storage_SHA_similarities(1,counter)=OEW_tentative; 
data_storage_SHA_similarities(2,counter)=OEW_similarities;

error=false;

% If the convergence criterion is not fulfilled, the error is set back to
% "true"
if OEW_similarities_deviation>0.005
       error=true;
       
       % Depending on the value, a different adjustment is made
       if OEW_similarities>OEW_tentative 
            MTOW_SHA_similarities=MTOW_SHA_similarities+1;
       else
            MTOW_SHA_similarities=MTOW_SHA_similarities-1;
       end
       
end % End of the conditional that checks the convergence

end % End of the loop

% Graphic plotting 
fig2=figure(2);
set(fig2,'Renderer', 'painters', 'Position', [400 400 500 350]);
hold on

plot([1:1:size(data_storage_SHA_similarities,2)], data_storage_SHA_similarities(1,:),'b')
hold on
plot([1:1:size(data_storage_SHA_similarities,2)], data_storage_SHA_similarities(2,:),'r')
legend('estimated','similarities')

% Axis format
set(gca,'TickLabelInterpreter','latex','fontsize',13)
ylabel('$OEW\,\left[\mathrm{kg}\right]$','interpreter','latex','FontSize',15)
xlabel('Number of iterations','interpreter','latex','FontSize',15)
xlim([0 1.5e7])

% Grid format
grid on
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.2;
grid minor

legend('Tentative','Similarities','Location','southeast','interpreter','latex','fontsize',13)