% Function to solve the Breguet's differential equation in terms of range
% using a 1st order Euler numerical method

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

function [X_sol,W_sol,frac,Efficiency] = Euler_range_function(c_t,M,C_D0,k,S,H,W_0,Range,N)

%% 1. Definition of Constants, Parameters and Variables

% 1.1. CONSTANTS   
Ru = 8.31432;           % Universal Constant for Ideal Gases    [J/mole*K]
% Earth
g = 9.80665;            % Acceleration at Earth's surface       [m/s^2]
T0 = 288.15;            % US Standard Sea Level Temperature     [K]
P0 = 101325;            % Pressure at Sea Level                 [Pa]
Mm = 28.9644*10^-3;     % Molecular Mass                        [kg*mole^-1]    
H_layer = 1e3*[0 11 20 32 47 52 61 69 79 90 100 110 117.776];   % Earth's atmospheric layers
lambda = 1e-3*[-6.5 0 1 2.8 0 -2 -4 -3 0 2 4.36 16.4596 0];     % Earth's atmospheric layers altitude thermal gradient [k/m]
gamma = 1.4;            % Earth's air specific heats relation   [adim]
R = Ru/Mm;              % Gas constant for Earth's air

%% 2. PHYSICAL DATA

% Atmospherical properties
[T_h,rho]=ISA_atmosphere(H);

%% Euler approximation

% 3.1 Numerical data
x_0 = 0;
x_f = Range;
h = (x_f-x_0)/N;

% Declaration of solution vectors
X_sol = linspace(x_0,x_f,N);
W_sol = zeros(1, length(X_sol));
Efficiency = zeros(1, length(X_sol));

% 3.2 Initial conditions
W_sol(1) = W_0;
X_sol(1) = x_0;

% 3.3 Euler E1
for i=1:1:length(X_sol)
    
    slope = -c_t*g/(M*sqrt(gamma*R*T_h))*(0.5*C_D0*rho*(M^2*gamma*R*T_h)*S+k*(2*W_sol(i)^2)/(rho*(M^2*gamma*R*T_h)*S));
    X_sol(i+1) = X_sol(i) + h;
    W_sol(i+1) = W_sol(i)+h*slope;    
    
    C_L=(2*W_sol(i))/(rho*(M^2*gamma*R*T_h)*S);
    C_D=C_D0+k*C_L^2;
    Efficiency(i)=C_L/C_D;
    
end % End of the Euler loop

% Fraction between final and initial mass
frac = W_sol(end)/W_0;

end

