function [C_D0,k] = Polar_parameters_function(W_0,S,Wing_span)

% Parasitic drag constant

a = -2.5229;
b = 1;
% Twin Engine Propellet Driven (1)
c_1 = 0.8635;
d_1 = 0.5632;
% Business Jet (2)
c_2 = 0.2263;
d_2 = 0.6977;
% Transport Jet (3)
c_3 = 0.0199;
d_3 = 0.7531;
% Imperial units
W_0_lb = W_0*2.20462; 
S_ft2 = S*(unitsratio('feet','meter'))^2; 

syms f1 f2 f3
f1 = solve (log10(f1) == a+b*(c_1+d_1*log10(W_0_lb)));
f2 = solve (log10(f2) == a+b*(c_2+d_2*log10(W_0_lb)));
f3 = solve (log10(f3) == a+b*(c_3+d_3*log10(W_0_lb)));

C_D0_i(1) = double(f1/S_ft2);                 
C_D0_i(2) = double(f2/S_ft2);
C_D0_i(3) = double(f3/S_ft2);

C_D0 = min(C_D0_i);

% Lift induced drag constant

A = Wing_span^2/S;  % Alargamiento 
e = 0.85;         % Factor de eficiencia
k = 1/(pi*A*e);