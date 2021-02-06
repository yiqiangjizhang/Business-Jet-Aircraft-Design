% Function to solve the Breguet's differential equation in terms of endurance
% using a 4th order Runge-Kutta numerical method

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

function [T_sol,W_sol,frac] = RK4_autonomy_function(c_t,M,C_D0,k,S,H,W_0,Autonomy)

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

% Atmospherical conditions
[T_h,rho]=ISA_atmosphere(H);

%% 3. RUNGE-KUTTA RK4

% 3.1 Numerical data
DeltaT = 10;

% Declaration of solution vectors
T_sol = 0:DeltaT:Autonomy;
W_sol = zeros(1, length(T_sol));

% 3.2 Initial conditions
W_sol(1) = W_0;
%[C_D0,k] = Polar_parameters_function(W_sol(1)/9.81,S,Wing_span);

% 3.3 Runge-Kutta RK4
for i = 1:length(T_sol)-1
    
    % Functions
    F1 = @(T, W) -c_t*(0.5*C_D0*rho*(M^2*gamma*R*T_h)*S+k*(2*W^2)/(rho*(M^2*gamma*R*T_h)*S));
    % Computation of coefficients sub 1
    i1 = F1(T_sol(i), W_sol(i));
    % Computation of coefficients sub 2
    i2 = F1(T_sol(i)+ DeltaT/2, W_sol(i)+i1*DeltaT/2);
    % Computation of coefficients sub 3
    i3 = F1(T_sol(i)+ DeltaT/2, W_sol(i)+i2*DeltaT/2);
    % Computation of coefficients sub 4
    i4 = F1(T_sol(i)+ DeltaT, W_sol(i)+i3*DeltaT);
    % Compute next step
    W_sol(i+1) = W_sol(i) + (DeltaT/6)*(i1 + 2*i2 + 2*i3 + i4);
    
end % End of the Runge-Kutta loop

% Fraction between final and initial mass
frac = W_sol(end)/W_0;

end






