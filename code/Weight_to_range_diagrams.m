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
%pax = 10 ; % Number of passengers
W_pax = 77 ; % Mass of a passenger [kg]
W_bagg = 20 ; % Mass of a passenger's baggage [kg]
rho_cargo = 160 ; % Cargo density [kg/m^3]
rho_bagg = 200 ;  % Baggage density [kg/m^3]
kbd = 0.85 ; % Hold occupation efficiency

% Performance data
% Expected total range

M_cruise = 0.8; % Cruise Mach number
Efficiency = 15 ; % Aerodynamic efficiency during the cruise
h_cruise = 12000 ; % Cruise height [m]
cj_cruise = 2*1.8 * 10^(-5); % Cruise specific consumption [kg/N/s]
cj_loiter = 2*1.1 * 10^(-5); % Loiter specific consumption [kg/N/s]
E_loiter = 30 * 60 ; % Loiter time [s]
%E_reserve= 30 * 60 ; % Reserved time to reach alternate aerodrome [s]
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


N=100;
M=500;
TOW=zeros(N,1);
OEW=zeros(N,1);
PL=zeros(N,1);
W_f_reserve_analytical_vec=zeros(N,1);
%MTOW_SHA_torenbeek=32626;
Range_vec=linspace(4500,5500,N);
total_range=Range_vec*unitsratio('meter','nauticalmile');  % [m]
R_cruise = 0.8*total_range; % Cruise range [m] 
PL_fullfuel=convmass(1700,'lbm','kg');


% pax=linspace(10,0,M);
pax=10;
% Baggage volume
V_bagg=k_b_baggage*b_f_expected^2*l_f_expected;     
% Maximum payload calculation
MPL = pax * (W_pax+W_bagg) + (kbd*V_bagg -((W_bagg*pax)/rho_bagg))*rho_cargo ;  % Computation fo Payload Weight
PL_vec=linspace(MPL,0,M);


E_reserve=linspace(30,5,M)*60; %[s]


converge_final=true;
for r=1:N
    if converge_final==false
        break;
    else
    R_cruise_2=R_cruise(r);
    for i=1:M
        %converge_anterior=converge;
        PL_2=PL_vec(i);
        E_reserve_2=E_reserve(i);
        [OEW_torenbeek,MTOW_SHA_torenbeek,W_f_analytical,W_f_trip,W_f_reserve_analytical]= Weight_calculation_ANALYTICAL_roskam(PL_2,R_cruise_2,E_reserve_2,b_f_expected,a_f_expected,l_f_expected,M_cruise,Efficiency,cj_cruise,cj_loiter,E_loiter,MTOW_estimated,engine_weight,T_cruise,R_g,lambda,g,Delta_W_e_regression);
        %if (i<=M && converge==false)
            if r==1
                OEW_bueno=OEW_torenbeek;
                OEW(r)=OEW_bueno;
                MTOW_MPL=MTOW_SHA_torenbeek;
                TOW(r)=MTOW_MPL;
                PL(r)=PL_2;
                W_f_reserve_analytical_vec(r)=W_f_reserve_analytical;
                break;
            elseif ( (abs(MTOW_SHA_torenbeek-MTOW_MPL)/MTOW_MPL)<0.0005 && (abs(OEW_torenbeek-OEW_bueno)/OEW_bueno)<0.0005 && PL_2>=PL_fullfuel) % MTOW condtion and maximum fuel weight payload condition
                TOW(r)=MTOW_SHA_torenbeek;
                OEW(r)=OEW_torenbeek;
                W_f_reserve_analytical_vec(r)=W_f_reserve_analytical;
                E_reserve_fullfuel=E_reserve_2;
                PL(r)=PL_2;
                converge=true;
                r_MTOW=r;
                i_Rmax=i;
                R_MTOW=Range_vec(r);
                W_f_analytical_max=W_f_analytical;
                W_f_reserve_fullfuel=W_f_reserve_analytical;
                break;
            else
                errorMTOW=abs(MTOW_SHA_torenbeek-MTOW_MPL)/MTOW_MPL;
                errorOEW=abs(OEW_torenbeek-OEW_bueno)/OEW_bueno;
                converge=false;

            end
        if (i==M && converge==false)
            converge_final=false;
            break;
        end
    end
    end
    
end


TOW_final=OEW_bueno+W_f_analytical_max;
M_ff_analytical_final=1-W_f_analytical_max/TOW_final ;



fraction1 = 0.990 ; % Mass fuel fraction engine start phase
fraction2 = 0.995 ; % Mass fuel fraction taxi phase
fraction3 = 0.995 ; % Mass fuel fraction take off phase
fraction4 = 0.980 ; % Mass fuel fraction Climb phase
fraction6_analytical = exp(-E_loiter/(Efficiency/(g*cj_loiter))); % Mass fuel fraction loiter phase analytical computation
fraction_reserve=exp(-E_reserve_fullfuel/(Efficiency/(g*cj_loiter))); % Mass fuel fraction for reaching alternate aerodrome
fraction7 = 0.990 ; % Mass fuel fraction descend phase
fraction8 = 0.992; % Mass fuel fraction shutdown phase


fraction5_analytical_final=M_ff_analytical_final/(fraction1*fraction2*fraction3*fraction4*fraction6_analytical*fraction7*fraction8*fraction_reserve);

ktf = Efficiency * M_cruise*sqrt(lambda*T_cruise*R_g)/(g*cj_cruise);

R_cruise_final=-log(fraction5_analytical_final)*ktf;
R_max=(R_cruise_final/0.8)*unitsratio('nauticalmile','meter');


Range_vec_1=[0 4500];
TOW_vec_1=[OEW_bueno+MPL+W_f_reserve_analytical_vec(1) MTOW_MPL];

Range_vec_2=[R_MTOW R_max];
TOW_vec_2=[TOW(r_MTOW) TOW_final];

PL_vec_1=[MPL MPL];
PL_vec_2=[PL(r_MTOW) 0];

W_f_reserve_analytical_vec_1=[W_f_reserve_analytical_vec(1) W_f_reserve_analytical_vec(1)];
W_f_reserve_analytical_vec_2=[W_f_reserve_analytical_vec(r_MTOW) W_f_reserve_analytical_vec(r_MTOW)];




%% DIAGRAMS

fig1=figure(1);
set(fig1,'Renderer', 'painters', 'Position', [400 400 800 450]);
hold on

plot(Range_vec_1,TOW_vec_1,'b','HandleVisibility','off');
hold on

plot(Range_vec(1:r_MTOW),TOW(1:r_MTOW),'b','DisplayName','$\mathit{TOW}$');
hold on

plot(Range_vec_2,TOW_vec_2,'b','HandleVisibility','off');
hold on

yline(OEW_bueno,'g','DisplayName','$\mathit{OEW}$');
hold on
xline(Range_vec(1),'--k','linewidth',1.2,'HandleVisibility','off');
hold on
xline(Range_vec(r_MTOW),'--k','linewidth',1.2,'HandleVisibility','off');
hold on
xline(Range_vec_2(2),'--k','linewidth',1.2,'HandleVisibility','off');

hold on
plot(Range_vec_1,OEW_bueno+PL_vec_1,'m','HandleVisibility','off');
hold on
plot(Range_vec(1:r_MTOW),OEW_bueno+PL(1:r_MTOW),'m','DisplayName','$\mathit{OEW}+PL$');
hold on
plot(Range_vec_2,OEW_bueno+PL_vec_2,'m','HandleVisibility','off');
hold on

plot(Range_vec_1,OEW_bueno+PL_vec_1+W_f_reserve_analytical_vec_1,'r','HandleVisibility','off');
hold on
plot(Range_vec(1:r_MTOW),OEW_bueno+PL(1:r_MTOW)+W_f_reserve_analytical_vec(1:r_MTOW),'r','DisplayName','$\mathit{OEW}+PL+RF$');
hold on
plot(Range_vec_2,OEW_bueno+PL_vec_2+W_f_reserve_analytical_vec_2,'r','HandleVisibility','off');
ylim([1E4 3.5e4])
xlim([3500 5000])


% Axis format
set(gca,'TickLabelInterpreter','latex','fontsize',15)
ylabel('Weight $\left[\mathrm{kg}\right]$','interpreter','latex','FontSize',16)
xlabel('Range $\left[\mathrm{NM}\right]$','interpreter','latex','FontSize',16)

% xticks([3500 4000]);
% xticklabels({'3500','4000'});

% Grid format
grid on
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.2;
grid minor


%
%text(4420,1.7E4,'$R_{\mathit{MPL}}=4500$','Color','k','Rotation',90,'interpreter','latex','Fontsize',16) ;
%text(4700,1.7E4,'$R_{MTOW}=4783$','Color','k','Rotation',90,'interpreter','latex','Fontsize',16) ; 
text(4460,1.7E4,'$R_{\mathit{MPL}}=4500$','Color','k','Rotation',90,'interpreter','latex','Fontsize',16) ; 
text(4740,1.7E4,'$R_{MTOW}=4783$','Color','k','Rotation',90,'interpreter','latex','Fontsize',16) ; 
text(4950,1.7E4,'$R_{\mathrm{max}}=4989$','Color','k','Rotation',90,'interpreter','latex','Fontsize',16) ; 

%text(4670,1.1E4,'$R_{MTOW}$','Color','k','interpreter','latex','Fontsize',14) ; 
%text(4900,1.1E4,'$R_{MAX}$','Color','k','interpreter','latex','Fontsize',14) ; 
legend('Location','Northeastoutside','Interpreter','latex');

xlim([3500 5000])

 %print(fig1,'diagrama_pesos_alcance1.png','-dpng','-r800'); 















fig2=figure(2);
set(fig2,'Renderer', 'painters', 'Position', [400 400 800 450]);
hold on

plot(Range_vec_1,TOW_vec_1,'b','HandleVisibility','off');
hold on

plot(Range_vec(1:r_MTOW),TOW(1:r_MTOW),'b','DisplayName','$TOW$');
hold on

plot(Range_vec_2,TOW_vec_2,'b','HandleVisibility','off');
hold on

yline(OEW_bueno,'g','DisplayName','$\mathit{OEW}$');
hold on
xline(Range_vec(1),'--k','linewidth',1.2,'HandleVisibility','off');
hold on
xline(Range_vec(r_MTOW),'--k','linewidth',1.2,'HandleVisibility','off');
hold on
xline(Range_vec_2(2),'--k','linewidth',1.2,'HandleVisibility','off');

hold on
plot(Range_vec_1,OEW_bueno+PL_vec_1,'m','HandleVisibility','off');
hold on
plot(Range_vec(1:r_MTOW),OEW_bueno+PL(1:r_MTOW),'m','DisplayName','$\mathit{OEW}+PL$');
hold on
plot(Range_vec_2,OEW_bueno+PL_vec_2,'m','HandleVisibility','off');
hold on

plot(Range_vec_1,OEW_bueno+PL_vec_1+W_f_reserve_analytical_vec_1,'r','HandleVisibility','off');
hold on
plot(Range_vec(1:r_MTOW),OEW_bueno+PL(1:r_MTOW)+W_f_reserve_analytical_vec(1:r_MTOW),'r','DisplayName','$\mathit{OEW}+PL+RF$');
hold on
plot(Range_vec_2,OEW_bueno+PL_vec_2+W_f_reserve_analytical_vec_2,'r','HandleVisibility','off');
ylim([0 35E3])
xlim([0 5000])

% Axis format
set(gca,'TickLabelInterpreter','latex','fontsize',15)
ylabel('Weight $\left[\mathrm{kg}\right]$','interpreter','latex','FontSize',16)
xlabel('Range $\left[\mathrm{NM}\right]$','interpreter','latex','FontSize',16)
 
% Grid format
grid on
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.2;
grid minor

text(4420,1.8E4,'$R_{\mathit{MPL}}=4500$','Color','k','Rotation',90,'interpreter','latex','Fontsize',15) ; 
text(4690,1.8E4,'$R_{MTOW}=4783$','Color','k','Rotation',90,'interpreter','latex','Fontsize',15) ; 
text(4905,1.8E4,'$R_{\mathrm{max}}=4989$','Color','k','Rotation',90,'interpreter','latex','Fontsize',15) ; 
legend('Location','Northwest','Interpreter','latex');
% 
% 
%print(fig2,'diagrama_pesos_alcance2.png','-dpng','-r800');












fig3=figure(3);
set(fig3,'Renderer', 'painters', 'Position', [400 400 800 500]);

hold on
xline(Range_vec(1),'--k','linewidth',1.2,'HandleVisibility','off');
hold on
xline(Range_vec(r_MTOW),'--k','linewidth',1.2,'HandleVisibility','off');
hold on
xline(Range_vec_2(2),'--k','linewidth',1.2,'HandleVisibility','off');

hold on
plot(Range_vec_1,PL_vec_1,'m','HandleVisibility','off');
hold on
plot(Range_vec(1:r_MTOW),PL(1:r_MTOW),'m','DisplayName','$PL$');
hold on
plot(Range_vec_2,PL_vec_2,'m','HandleVisibility','off');
hold on

ylim([0 1.4E3])
xlim([0 5000])

% Axis format
set(gca,'TickLabelInterpreter','latex','fontsize',14)
ylabel('Payload $\left[\mathrm{kg}\right]$','interpreter','latex','FontSize',15)
xlabel('Range $\left[\mathrm{NM}\right]$','interpreter','latex','FontSize',15)
 
% Grid format
grid on
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.2;
grid minor

text(4410,0.01E4,'$R_{\mathit{MPL}}=4500$','Color','k','Rotation',90,'interpreter','latex','Fontsize',16) ; 
text(4690,0.01E4,'$R_{MTOW}=4783$','Color','k','Rotation',90,'interpreter','latex','Fontsize',16) ; 
text(4900,0.09E4,'$R_{\mathrm{max}}=4989$','Color','k','Rotation',90,'interpreter','latex','Fontsize',16) ; 
%legend('Location','Southwest','Interpreter','latex');

%print(fig3,'payload_to_range_1.png','-dpng','-r800');


% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
fig4=figure(4);
set(fig4,'Renderer', 'painters', 'Position', [400 400 800 500]);

hold on
xline(Range_vec(1),'--k','linewidth',1.2,'HandleVisibility','off');
hold on
xline(Range_vec(r_MTOW),'--k','linewidth',1.2,'HandleVisibility','off');
hold on
xline(Range_vec_2(2),'--k','linewidth',1.2,'HandleVisibility','off');

hold on
plot(Range_vec_1,PL_vec_1,'m','HandleVisibility','off');
hold on
plot(Range_vec(1:r_MTOW),PL(1:r_MTOW),'m','DisplayName','PL');
hold on
plot(Range_vec_2,PL_vec_2,'m','HandleVisibility','off');
hold on

ylim([0 1.4E3])
xlim([3500 5000])

% Axis format
set(gca,'TickLabelInterpreter','latex','fontsize',14)
ylabel('Payload $\left[\mathrm{kg}\right]$','interpreter','latex','FontSize',15)
xlabel('Range $\left[\mathrm{NM}\right]$','interpreter','latex','FontSize',15)
 
% Grid format
grid on
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.2;
grid minor

text(4470,0.01E4,'$R_{\mathit{MPL}}=4500$','Color','k','Rotation',90,'interpreter','latex','Fontsize',16) ; 
text(4750,0.01E4,'$R_{MTOW}=4783$','Color','k','Rotation',90,'interpreter','latex','Fontsize',16) ; 
text(4955,0.07E4,'$R_{\mathrm{max}}=4989$','Color','k','Rotation',90,'interpreter','latex','Fontsize',16) ; 

%print(fig4,'payload_to_range_2.png','-dpng','-r800');

% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
