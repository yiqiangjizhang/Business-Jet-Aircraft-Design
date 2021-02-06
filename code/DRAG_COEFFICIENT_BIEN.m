% Code that computes the drag polars using the component buildup method

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

% Select the desired drag contribution
% Type 1: Cruise flight considering transonic effects
% Type 2: Take-off flaps effect
% Type 3: Landing flaps effect
% Type 4: Landing gear wheels effect
type_of_drag=1;

% Flight conditions

% Cruise Mach number
M_cruise = 0.8; 
% Cruise height
h_cruise=12000; % [m]
% ISA atmosphere conditons
[a,rho,mu]=ISA_atmosphere_drag(h_cruise); % [m/s, kg/m^3, kg/m*s]

% Fuselage dimensions
fuselage_length=21.4; % [m]
fuselage_diameter=2; % [m]
slenderness_fuselage=fuselage_length/fuselage_diameter;
interference_factor_fuselage=1;

% Engine nacelle dimensions
nacelle_length=128*unitsratio('meter','inch'); % [m]
nacelle_diameter=52*unitsratio('meter','inch'); % [m] 
slenderness_nacelle=nacelle_length/nacelle_diameter;
interference_factor_nacelle=1.5;

% Main wing dimensions
main_wing_surface=70.2; % [m^2]
mean_aerodynamic_chord_main_wing=3.19; % [m]
lambda_sweep_main_wing=deg2rad(27); % [rad]
wingspan_main_wing=23; % [m]
aspect_ratio_main_wing=wingspan_main_wing^2/main_wing_surface;
interference_factor_main_wing=1;
% SC 0402 airfoil
xc_main_wing=0.34;
tc_main_wing=0.02;
% Reference surface assignation
Reference_surface=main_wing_surface;

% Horizontal stabilizer
mean_aerodynamic_chord_horizontal_stabilizer=1.9879; % [m]
lambda_sweep_horizontal_stabilizer=deg2rad(27); % [rad]
wingspan_horizontal_stabilizer=7.4547; % [m]
interference_factor_horizontal_stabilizer=1.05;
% NACA 0009 airfoil
xc_horizontal_stabilizer=0.309;
tc_horizontal_stabilizer=0.09;

% Vertical stabilizer
mean_aerodynamic_chord_vertical_stabilizer=3.5464; % [m]
lambda_sweep_vertical_stabilizer=deg2rad(27); % [rad]
wingspan_vertical_stabilizer=3.5464; % [m]
interference_factor_vertical_stabilizer=1.05;
% NACA 0009 airfoil
xc_vertical_stabilizer=0.309;
tc_vertical_stabilizer=0.09;

% Flaps
MAC_ratio_flap=0.3;
A=0.0016;
B=1.5;
deflection_angle_take_off=25; % [degrees]
deflection_angle_landing=37; % [degrees]

% Landing gear wheels
% Main gear
number_wheels_main_gear=4;
diameter_wheels_main_gear=0.8636; % [m]
width_wheels_main_gear=0.235; % [m]
surface_wheels_main_gear=diameter_wheels_main_gear*width_wheels_main_gear; % [m^2]
% Nose gear
number_wheels_nose_gear=2;
diameter_wheels_nose_gear=0.5334; % [m]
width_wheels_nose_gear=0.1842; % [m]
surface_wheels_nose_gear=diameter_wheels_nose_gear*width_wheels_nose_gear; % [m^2]

%% INDUCED DRAG

% Oswald parameter (in this case, phi=e)
Oswald_parameter=4.61*(1-0.045*aspect_ratio_main_wing^(0.68))*(cos(lambda_sweep_main_wing))^0.15-3.1;

% Induced drag constant k
k_induced=1/(pi*aspect_ratio_main_wing*Oswald_parameter);

%% PARASITIC DRAG

% Reynolds function definition
Reynolds_function=@(characteristic_length) rho*M_cruise*a*characteristic_length/mu;
% Flat-plate skin-friction coefficient definition
skin_friction_coefficient=@(Reynolds_number) 0.455/((log10(Reynolds_number) )^2.58*(1+0.144*M_cruise^2)^0.65);

% WING AND TAIL

% Form factor function definition
FF_wingtail_function=@(xc,tc,Lambda_wing) 1.1*(1+(0.6/xc)*tc+100*tc^4)*(1.34*(M_cruise^0.18)*(cos(Lambda_wing))^0.28);
% Wet surface function definition
S_wet_wingtail_function=@(tc,b,MAC) 2*(1+0.5*tc)*b*MAC;

% MAIN WING
% Form factor
FF_main_wing=FF_wingtail_function(xc_main_wing,tc_main_wing,lambda_sweep_main_wing);
% Reynolds number
Reynolds_main_wing=Reynolds_function(mean_aerodynamic_chord_main_wing);
% Flat-plate skin-friction coefficient
skin_friction_coefficient_main_wing=skin_friction_coefficient(Reynolds_main_wing);
% Wet surface
S_wet_main_wing=S_wet_wingtail_function(tc_main_wing,wingspan_main_wing,mean_aerodynamic_chord_main_wing);
% Drag coefficient
drag_coefficient_main_wing=FF_main_wing*skin_friction_coefficient_main_wing*interference_factor_main_wing*S_wet_main_wing/Reference_surface;

% HORIZONTAL STABILIZER
% Form factor
FF_horizontal_stabilizer=FF_wingtail_function(xc_horizontal_stabilizer,tc_horizontal_stabilizer,lambda_sweep_horizontal_stabilizer);
% Reynolds number
Reynolds_horizontal_stabilizer=Reynolds_function(mean_aerodynamic_chord_horizontal_stabilizer);
% Flat-plate skin-friction coefficient
skin_friction_coefficient_horizontal_stabilizer=skin_friction_coefficient(Reynolds_horizontal_stabilizer);
% Wet surface
S_wet_horizontal_stabilizer=S_wet_wingtail_function(tc_horizontal_stabilizer,wingspan_horizontal_stabilizer,mean_aerodynamic_chord_horizontal_stabilizer);
% Drag coefficient
drag_coefficient_horizontal_stabilizer=FF_horizontal_stabilizer*skin_friction_coefficient_horizontal_stabilizer*interference_factor_horizontal_stabilizer*S_wet_horizontal_stabilizer/Reference_surface;

% VERTICAL STABILIZER
% Form factor
FF_vertical_stabilizer=FF_wingtail_function(xc_vertical_stabilizer,tc_vertical_stabilizer,lambda_sweep_vertical_stabilizer);
% Reynolds number
Reynolds_vertical_stabilizer=Reynolds_function(mean_aerodynamic_chord_vertical_stabilizer);
% Flat-plate skin-friction coefficient
skin_friction_coefficient_vertical_stabilizer=skin_friction_coefficient(Reynolds_vertical_stabilizer);
% Wet surface
S_wet_vertical_stabilizer=S_wet_wingtail_function(tc_vertical_stabilizer,wingspan_vertical_stabilizer,mean_aerodynamic_chord_vertical_stabilizer);
% Drag coefficient
drag_coefficient_vertical_stabilizer=FF_vertical_stabilizer*skin_friction_coefficient_vertical_stabilizer*interference_factor_vertical_stabilizer*S_wet_vertical_stabilizer/Reference_surface;


% FUSELAGE
% Form factor
FF_fuselage=0.9+5/((slenderness_fuselage)^1.5)+slenderness_fuselage/400;
% Reynolds number
Reynolds_fuselage=Reynolds_function(fuselage_length);
% Flat-plate skin-friction coefficient
skin_friction_coefficient_fuselage=skin_friction_coefficient(Reynolds_fuselage);
% Wet surface
S_wet_fuselage=pi*fuselage_length*fuselage_diameter;
% Drag coefficient
drag_coefficient_fuselage=FF_fuselage*skin_friction_coefficient_fuselage*interference_factor_fuselage*S_wet_fuselage/Reference_surface;

% ENGINE NACELLES
% Form factor
FF_nacelle=1+0.35/slenderness_nacelle;
% Reynolds number
Reynolds_nacelle=Reynolds_function(fuselage_length);
% Flat-plate skin-friction coefficient
skin_friction_coefficient_nacelle=skin_friction_coefficient(Reynolds_nacelle);
% Wet surface
S_wet_nacelle=pi*nacelle_diameter*nacelle_length;
% Drag coefficient
% Multiply times 2 because there two engines are mounted
drag_coefficient_nacelles=2*FF_nacelle*skin_friction_coefficient_nacelle*interference_factor_nacelle*S_wet_nacelle/Reference_surface;

% Final parasitic drag coefficient considering all the essential contributions
final_drag_coefficient=drag_coefficient_fuselage+drag_coefficient_main_wing+drag_coefficient_nacelles+drag_coefficient_vertical_stabilizer+drag_coefficient_horizontal_stabilizer;
% Final parasitic drag coefficient considering the leakages
final_drag_coefficient_with_leakages=1.02*final_drag_coefficient;
% Final parasitic drag coefficient considring the transsonic effects
final_drag_coefficient_with_transonic=1.4*final_drag_coefficient_with_leakages;

% FLAPS
% Drag coefficient function definition
drag_coefficient_flaps_function=@(deflection_angle) A*MAC_ratio_flap*(deflection_angle^B);
% Take-off drag coefficient obtention
drag_coefficient_flaps_take_off=drag_coefficient_flaps_function(deflection_angle_take_off);
% Landing drag coefficient obtention
drag_coefficient_flaps_landing=drag_coefficient_flaps_function(deflection_angle_landing);

% WHEELS
% Direct computation
drag_coefficient_wheels=0.3/Reference_surface*(number_wheels_main_gear*surface_wheels_main_gear+number_wheels_nose_gear*surface_wheels_nose_gear);

% Lift coefficient range used for the graphs
C_L_graph = [-2:0.01:3];

if type_of_drag==1 % Type 1: Cruise flight considering transonic effects

% Polar computed using Roskam's method    
C_D_cruise_Roskam = 0.020+0.0438*C_L_graph.^2;
% Polar computed thorugh component buildup method
C_D_cruise_buildup= final_drag_coefficient_with_transonic+k_induced*C_L_graph.^2;
    
% Figure size definition    
fig1=figure(1);
set(fig1,'Renderer', 'painters', 'Position', [400 400 500 350]);
hold on

% Polar plotting
plot(C_L_graph,C_D_cruise_buildup,'b','Linewidth',1)
plot(C_L_graph,C_D_cruise_Roskam,'r','Linewidth',1)

% Grid format
grid on;
grid minor;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.2;

% Axis format
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlabel('Lift coefficient $C_L$','interpreter','latex','FontSize',14)
ylabel('Drag coefficient $C_D$','interpreter','latex','FontSize',14)
legend('Buildup method','Roskam''s method','location','northwest','interpreter','latex');

% Printing command
%print(gcf,'final_cruise_drag.png','-dpng','-r700');
        
elseif type_of_drag==2 % Type 2: Take-off flaps effect
    
% Polar computed using Roskam's method        
C_D_take_off_Roskam = 0.035+0.0448*C_L_graph.^2;
% Polar computed thorugh component buildup method
C_D_take_off_buildup= final_drag_coefficient_with_leakages+drag_coefficient_flaps_take_off+k_induced*C_L_graph.^2;
    
% Figure size definition    
fig1=figure(1);
set(fig1,'Renderer', 'painters', 'Position', [400 400 500 350]);
hold on

% Polar plotting
plot(C_L_graph,C_D_take_off_buildup,'b','Linewidth',1)
plot(C_L_graph,C_D_take_off_Roskam,'r','Linewidth',1)

% Grid format
grid on;
grid minor;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.2;

% Axis format
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlabel('Lift coefficient $C_L$','interpreter','latex','FontSize',14)
ylabel('Drag coefficient $C_D$','interpreter','latex','FontSize',14)
legend('Buildup method','Roskam''s method','location','northwest','interpreter','latex');
    
% Printing command
%print(gcf,'final_takeoff_drag.png','-dpng','-r700');    
    
elseif type_of_drag==3 % Type 3: Landing flaps effect

% Polar computed using Roskam's method    
C_D_landing_Roskam = 0.085+0.0517*C_L_graph.^2;
% Polar computed thorugh component buildup method
C_D_landing_buildup= final_drag_coefficient_with_leakages+drag_coefficient_flaps_landing+k_induced*C_L_graph.^2;
    
% Figure size definition    
fig1=figure(1);
set(fig1,'Renderer', 'painters', 'Position', [400 400 500 350]);
hold on

% Polar plotting
plot(C_L_graph,C_D_landing_buildup,'b','Linewidth',1)
plot(C_L_graph,C_D_landing_Roskam,'r','Linewidth',1)

% Grid format
grid on;
grid minor;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.2;

% Axis format
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlabel('Lift coefficient $C_L$','interpreter','latex','FontSize',14)
ylabel('Drag coefficient $C_D$','interpreter','latex','FontSize',14)
legend('Buildup method','Roskam''s method','location','northwest','interpreter','latex');
    
% Printing command
%print(gcf,'final_landing_drag.png','-dpng','-r700');    
    
else % Type 4: Landing gear wheels effect

% Polar computed using Roskam's method    
C_D_landing_Roskam = 0.040+0.0438*C_L_graph.^2;
% Polar computed thorugh component buildup method
C_D_landing_buildup= final_drag_coefficient_with_leakages+drag_coefficient_wheels+k_induced*C_L_graph.^2;
    
% Figure size definition    
fig1=figure(1);
set(fig1,'Renderer', 'painters', 'Position', [400 400 500 350]);
hold on

% Polar plotting
plot(C_L_graph,C_D_landing_buildup,'b','Linewidth',1)
plot(C_L_graph,C_D_landing_Roskam,'r','Linewidth',1)

% Grid format
grid on;
grid minor;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.2;

% Axis format
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlabel('Lift coefficient $C_L$','interpreter','latex','FontSize',14)
ylabel('Drag coefficient $C_D$','interpreter','latex','FontSize',14)
legend('Buildup method','Roskam''s method','location','northwest','interpreter','latex');
    
% Printing command
%print(gcf,'final_wheel_drag.png','-dpng','-r700');    

end