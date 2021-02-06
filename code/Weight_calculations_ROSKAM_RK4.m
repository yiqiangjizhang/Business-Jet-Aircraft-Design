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
M_loiter = 0.6; % Cruise Mach number
R_cruise = 0.8*total_range; % Cruise range [m] 

% Cruise height vector. All the heights will be tackled equitatively
h_cruise_vector=[12000]; % [m]
% Length of each phase of the cruise
R_stage=R_cruise/length(h_cruise_vector); % [m]

h_loiter = 10000 ; % Loiter height [m]
h_reserve = 450; % Height for the reserve flight, stablished by ICAO [m]
cj_cruise = 2*1.8 * 10^(-5); % Cruise specific consumption [kg/N/s]
cj_loiter = 2*1.1 * 10^(-5); % Loiter specific consumption [kg/N/s]
E_loiter = 20 * 60 ; % Loiter time [s]
Autonomy= 30 * 60 ; % Reserved time to reach alternate aerodrome [s]
MTOW_estimated = 30100 ; % MTOW estimation [kg]
engine_weight=634.225; % Average engine weight [kg]

% Termo-physical properties
R_g = 2.870528738362446e+02 ; % Constant for air
lambda = 1.4 ; % Heat relation for air
g = 9.80665 ; % Gravity acceleration [m/s^2]

S = 51.1; % Wing surface [m^2]
Wing_span=20.9; % Wing span [m]

% Data regression loading
load('k_b_baggage.mat'); % Constant for the baggage volume
load('alpha.mat'); % Constant for the similarities criterion
load('Delta_W_e_regression.mat'); % Term for the Torenbeek criterion

% Number of points to numerically study the cruise
N=1e6;
% Number of points for each phase of the cruise
M=round(N/length(h_cruise_vector));

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
fraction7 = 0.990 ; % Mass fuel fraction descend phase
fraction8 = 0.992; % Mass fuel fraction shutdown phase

%% TORENBEEK CRITERION LOOP FOLLOWING ROSKAM'S METHOD


% Delta W_e abscissa calculation
x_fuselage = l_f_expected * 0.5 * (a_f_expected+b_f_expected);
% Delta_e computation
Delta_e = 10^Delta_W_e_regression(2) * x_fuselage^(Delta_W_e_regression(1));

% Initial error definition
error=true;
% Assignation of the estimated MTOW
MTOW_RK4_torenbeek=MTOW_estimated;
% Counter initialization
counter=0;
% Data storage matrix initialization
data_storage_RK4_torenbeek=[];

while error==true % Loop that is executed until convergence is achieved
    
counter=counter+1; % Counter addition

% Drag polar coefficients computation.
% It is only performed once, as it only depends on the MTOW
[C_D0,k] = Polar_parameters_function(MTOW_RK4_torenbeek,S,Wing_span);

% Mass at the beginning of the cruise stage
Mass_begin_cruise=fraction1*fraction2*fraction3*fraction4*MTOW_RK4_torenbeek;

% Fraction 5 initialization
fraction5_prod=1;

for stage=1:length(h_cruise_vector) % Loop that goes all over the phases of the cruise

% Cruise height assignation
h_cruise=h_cruise_vector(stage);
 
% Runge-Kutta resolution of this phase of the cruise
[~,~,fraction5(stage),Efficiency(stage,:)] = RK4_range_function(cj_cruise,M_cruise,C_D0,k,S,h_cruise,Mass_begin_cruise*9.81,R_stage,M);

% Mass at the end of the current stage =
% mass at the beginning of the next stage
Mass_begin_cruise=Mass_begin_cruise*fraction5(stage);

% Fraction 5 productory for phases
fraction5_prod=fraction5_prod*fraction5(stage);

end % End of the cruise loop

% The mass at the beginning of the loiter is the mass at the end of the
% cruise stage (i.e. the last cruise phase)
Mass_begin_loiter=Mass_begin_cruise;

% Runge-Kutta resolution of the loiter
[~,~,fraction6] = RK4_autonomy_function(cj_loiter,M_loiter,C_D0,k,S,h_loiter,Mass_begin_loiter*9.81,Autonomy);

% Mass at the beginning of the reserve stage
Mass_begin_reserve=Mass_begin_loiter*fraction6;

% Runge-Kutta resolution of the reserve stage
[~,~,fraction_reserve] = RK4_autonomy_function(cj_loiter,M_loiter,C_D0,k,S,h_reserve,Mass_begin_reserve*9.81,Autonomy);

% Mission fuel fraction
M_ff_RK4=fraction1*fraction2*fraction3*fraction4*fraction5_prod*fraction6*fraction_reserve*fraction7*fraction8

% Total fuel weight computation
W_f_RK4 = (1-M_ff_RK4)*MTOW_RK4_torenbeek; 

% Operating empty weight computation (tentative)
OEW_tentative=MTOW_RK4_torenbeek-W_f_RK4-MPL;

% Allowable operating empty weight computation (according to Torenbeek criterion)
OEW_torenbeek = 0.2*MTOW_RK4_torenbeek + 2*engine_weight+ Delta_e+500;

% Relative deviation computation
OEW_torenbeek_deviation=abs(OEW_torenbeek-OEW_tentative)/OEW_torenbeek;

% Data storage
data_storage_RK4_torenbeek(1,counter)=OEW_tentative; % Tentative result
data_storage_RK4_torenbeek(2,counter)=OEW_torenbeek; % Torenbeek result

error=false;

% If the convergence criterion is not fulfilled, the error is set back to
% "true"
if OEW_torenbeek_deviation>0.0005
       error=true;
       
       % Depending on the value, a different adjustment is made
       if OEW_torenbeek>OEW_tentative 
            MTOW_RK4_torenbeek=MTOW_RK4_torenbeek+1;
       else
            MTOW_RK4_torenbeek=MTOW_RK4_torenbeek-1;
       end
       
end % End of the conditional that checks the convergence

end % End of the loop

% Graphic plotting 
fig1=figure(1);
set(fig1,'Renderer', 'painters', 'Position', [400 400 500 350]);
hold on

plot([1:1:size(data_storage_RK4_torenbeek,2)], data_storage_RK4_torenbeek(1,:),'b')

plot([1:1:size(data_storage_RK4_torenbeek,2)], data_storage_RK4_torenbeek(2,:),'r')

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

legend('Tentative','Torenbeek','Location','northeast','interpreter','latex','fontsize',13)

%% SIMILARITIES CRITERION LOOP FOLLOWING ROSKAM'S METHOD

% Initial error definition
error=true;
% Assignation of the estimated MTOW
MTOW_RK4_similarities=MTOW_estimated;
% Counter initialization
counter=0;
% Data storage matrix initialization
data_storage_RK4_similarities=[];

while error==true % Loop that is executed until convergence is achieved
    
counter=counter+1; % Counter addition

% Drag polar coefficients computation.
% It is only performed once, as it only depends on the MTOW
[C_D0,k] = Polar_parameters_function(MTOW_RK4_similarities,S,Wing_span);

% Mass at the beginning of the cruise stage
Mass_begin_cruise=fraction1*fraction2*fraction3*fraction4*MTOW_RK4_similarities;

% Fraction 5 initialization
fraction5_prod=1;

for stage=1:length(h_cruise_vector) % Loop that goes all over the phases of the cruise

% Cruise height assignation
h_cruise=h_cruise_vector(stage);
 
% Runge-Kutta resolution of this phase of the cruise
[~,~,fraction5(stage),Efficiency(stage,:)] = RK4_range_function(cj_cruise,M_cruise,C_D0,k,S,h_cruise,Mass_begin_cruise*9.81,R_stage,M);

% Mass at the end of the current stage =
% mass at the beginning of the next stage
Mass_begin_cruise=Mass_begin_cruise*fraction5(stage);

% Fraction 5 productory for phases
fraction5_prod=fraction5_prod*fraction5(stage);

end % End of the cruise loop

% The mass at the beginning of the loiter is the mass at the end of the
% cruise stage (i.e. the last cruise phase)
Mass_begin_loiter=Mass_begin_cruise;

% Runge-Kutta resolution of the loiter
[~,~,fraction6] = RK4_autonomy_function(cj_loiter,M_loiter,C_D0,k,S,h_loiter,Mass_begin_loiter*9.81,Autonomy);

% Mass at the beginning of the reserve stage
Mass_begin_reserve=Mass_begin_loiter*fraction6;

% Runge-Kutta resolution of the reserve stage
[~,~,fraction_reserve] = RK4_autonomy_function(cj_loiter,M_loiter,C_D0,k,S,h_reserve,Mass_begin_reserve*9.81,Autonomy);

% Mission fuel fraction
M_ff_RK4=fraction1*fraction2*fraction3*fraction4*fraction5_prod*fraction6*fraction_reserve*fraction7*fraction8

% Total fuel weight computation
W_f_RK4 = (1-M_ff_RK4)*MTOW_RK4_similarities; 

% Operating empty weight computation (tentative)
OEW_tentative=MTOW_RK4_similarities-W_f_RK4-MPL;

% Allowable operating empty weight computation (according to similarities criterion)
OEW_similarities= alpha*MTOW_RK4_similarities;

% Relative deviation computation
OEW_similarities_deviation=abs(OEW_similarities-OEW_tentative)/OEW_similarities;

% Data storage
data_storage_RK4_similarities(1,counter)=OEW_tentative; % Tentative result
data_storage_RK4_similarities(2,counter)=OEW_similarities; % Similarities result

error=false;

% If the convergence criterion is not fulfilled, the error is set back to
% "true"
if OEW_similarities_deviation>0.0005
       error=true;
       
       % Depending on the value, a different adjustment is made
       if OEW_similarities>OEW_tentative 
            MTOW_RK4_similarities=MTOW_RK4_similarities+10;
       else
            MTOW_RK4_similarities=MTOW_RK4_similarities-10;
       end
       
end % End of the conditional that checks the convergence

end % End of the convergence loop

% Graphic plotting 
fig2=figure(2);
set(fig2,'Renderer', 'painters', 'Position', [400 400 500 350]);
hold on

plot([1:1:size(data_storage_RK4_similarities,2)], data_storage_RK4_similarities(1,:),'b')
hold on
plot([1:1:size(data_storage_RK4_similarities,2)], data_storage_RK4_similarities(2,:),'r')
legend('estimated','similarities')

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

legend('Tentative','Similarities','Location','east','interpreter','latex','fontsize',13)