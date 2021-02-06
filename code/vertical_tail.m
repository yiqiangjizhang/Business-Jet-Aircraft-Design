%% Data entry

S=70.2; %wing reference area
b=23; %wingspan of the main wing
l_optim=9.5; %tail moment arm
V_v=0.074; %vertical tail volume coefficient
AR_v=1; %vertical tail aspect ratio
TR_v=0.8; %vertical tail taper ratio
i_v=0; %vertical tail incidence angle
sweep_angle_v=deg2rad(27); %vertical tail sweep angle
dihedral_angle_v=0; %vertical tail dihedral angle
area_f=pi; %mean fuselage area
%vertical tail profile NACA 0009
Cl_alpha_v=0.0968; %vertical tail profile Cl_alpha
Z_w=0.9; %distance wingroot to middle fuselage
eta_v=0.96; %vertical tail dynamic pressure ratio
K_f1=0.75; %fuselage contribution to C_n_beta

%% Calculation of the vertical tail planform area

S_v=b*S*V_v/l_optim; %vertical tail planform area

%% Calculation of the vertical tail parameters

b_v=sqrt(AR_v*S_v); %vertical tail wingspan
C_mean=b_v/AR_v; %MAC of the vertical tail
Cv_root=C_mean*(3*(1+TR_v))/(2*(1+TR_v+TR_v^2)); %root chord of the vertical tail
Cv_tip=Cv_root*TR_v; %tip chord of the vertical tail
CL_alpha_v=Cl_alpha_v/(1+(Cl_alpha_v)/(pi*AR_v)); %vertical tail CL_alpha

%% Lateral static stability criteria

d=sqrt(area_f/0.7854);
der_sigma=-0.276+3.06*S_v/(S*(1+cos(deg2rad(27))))+0.4*Z_w/d+0.009*AR_w; %sidewash
Cn_beta_v=K_f1*CL_alpha_v*(1-der_sigma)*eta_v*l_optim*S_v/(b*S); %vertical tail C_n_beta