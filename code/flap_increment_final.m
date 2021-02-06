%% Eta regression for plain flap, curved part
clc;
clear;
% Plot regression
eta=[0.8,0.8,0.74,0.69,0.52,0.45,  0.4,0.37];
delta=[5,10,17,20,30,40,50,60];
delta_long=5:2:60;
eta_coeff = polyfit(delta,eta,4);
figure(1)
plot(delta_long,polyval(eta_coeff,delta_long),delta,eta,'o'); 
legend('regression','data')

%% Calculation of delta CL
% Aircraft data
Cl_alpha=6.6463;
cf_c=0.3;
theta_f=acos(2*cf_c-1);
Tau=1-(theta_f-sin(theta_f))/pi;
sweep=deg2rad(27);
Sw=70.2;
b=23; %Wing span [m]
C_m=3.19; %Wing Mean chord [m]
Swf=49.15;% Flap wing surface [m]

% Lift increment calculation
delta_f=0:0.5:60;
for i=1:length(delta_f)
    eta_function(i)=polyval(eta_coeff,delta_f(i));
    if delta_f(i)<10
        delta_Cl0(i)=Cl_alpha*Tau*0.8*deg2rad(delta_f(i));
    else 
        delta_Cl0(i)=Cl_alpha*Tau*eta_function(i)*deg2rad(delta_f(i));
    end
    delta_Clmax(i)=delta_Cl0(i);
    delta_CLmax(i)=0.92*delta_Clmax(i)*Swf/Sw*cos(sweep); 
end

%Plot of increment of maximum lift coefficient
figure(2)
plot(delta_f,delta_CLmax)
xlabel('$$Flap\ deflection\ \delta_{f} [^o]$$','Interpreter','latex','FontSize',12) ;
ylabel('$$\Delta C_{L,MAX}$$','Interpreter','latex','FontSize',12) ; 
set(gca,'TickLabelInterpreter','latex')
title('Flap effect on maximum lift coefficient','Interpreter','latex','FontSize',15);
grid on
grid minor

%% Calculation of maximum delta CL for delta Cl=1.3
pos=find(delta_Clmax>=1.292 & delta_Clmax<=1.302);
fprintf('\n Angle %d',delta_f(pos));
fprintf('\n Increment max Cl %d',delta_Clmax(pos));
fprintf('\n Increment max CL %d',delta_CLmax(pos));