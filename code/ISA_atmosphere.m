function [T,rho]=ISA_atmosphere(Z)

% Physical constants definition
R_asterisc=8.31432; % [N m/(mol K)]
M_0=28.964420*1e-3; % [kg/mol] 
N_A=6.02257e26*1e-3; % [1/mol]
g_0_geo=9.80665; % [m^2/(s^2 m']
g_0=9.80665; % [m/s^2]
Gamma=g_0_geo/g_0; % [m'/m]
r_0=6.356766e6; % [m]

% Geopotential height vector
H_b=[0 11 20 32 47 51 71 80]*1e3; % [m']
% Temperature gradient vector
L_Mb=[-6.5 0.0 1.0 2.8 0.0 -2.8 -2.0]*1e-3; % [K/m']
beta_viscosa=1.458*10^(-6); % [kg /(s m K^1/2]
S=110.4; %[K] Sutherland's constant
gamma_aire=1.4;

% Computation of the base pressures (P_b) and the base temperatures at a
% mollecular scale (T_b) of each region of the atmosphere.
% Data obtained from ISO 2533:1975 (ISA atmosphere)
T_Mb=[288.15 zeros(1,length(L_Mb)-1)]; % [K]
P_b=[101325.0 zeros(1,length(L_Mb)-1)]; % [Pa]
for i=1:(length(L_Mb)-1)
    T_Mb(i+1)=T_Mb(i)+L_Mb(i)*(H_b(i+1)-H_b(i));
    if L_Mb(i)==0
        P_b(i+1)=P_b(i)*exp((-g_0_geo*M_0*(H_b(i+1)-H_b(i)))/(R_asterisc*T_Mb(i)));
    else
        P_b(i+1)=P_b(i)*(T_Mb(i)/(T_Mb(i)+L_Mb(i)*(H_b(i+1)-H_b(i))))^(g_0_geo*M_0/(R_asterisc*L_Mb(i)));
    end
end

% Geopotential altitude computation
H=Z*Gamma*r_0/(r_0+Z);

% Temperature and pressure computation according to the geopotential
% altitude according to ISO 2533:1975
for b=1:(length(L_Mb))
    if H>=H_b(b) && H<H_b(b+1) % Es distingeix el tram de l'atmosfera
       T_M=T_Mb(b)+L_Mb(b)*(H-H_b(b));
       if L_Mb(b)==0 % Si existeix gradient de temperatura
           P=P_b(b)*exp((-g_0_geo*M_0*(H-H_b(b)))/(R_asterisc*T_Mb(b)));
       else % Si no existeix gradient de temperatura
           P=P_b(b)*(T_Mb(b)/(T_Mb(b)+L_Mb(b)*(H-H_b(b))))^(g_0_geo*M_0/(R_asterisc*L_Mb(b)));
       end
       break
    end
end

% Gravity acceleration computation according to Z
g=g_0*(r_0./(r_0+Z)).^2;

% Sound speed calculation from T_M
a=sqrt(gamma_aire*R_asterisc*T_M/M_0);

% Other termophysical properties
T=T_M;
rho=(P*M_0)./(T*R_asterisc); 
mu=(beta_viscosa*T.^(1.5))./(S+T); 

end