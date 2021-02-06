%% 1. Data entry.

m_chord=3.19; %wing MAC
S=70.2; %wing reference area
b=23; %wingspan
L_f=21.4; %fuselage length
D_f=2.0; %maximum fuselage diameter
Cl_w=0.33; %CL cruise
alpha_t=0; %angle incidencia fuselatge
AR_w=7.5; %aspect ratio wing
TR_w=0.35; %tapper ratio wing
sweep_angle=deg2rad(27); %sweep angle
twist_angle=deg2rad(-2); %twist angle
Cm_o_w=-0.044; %wing cm_o
h_o=0.25; %wing ca location relative to wing MAC
dist_Cg_Ca_w=-0.027; %cg to wing ca distance
cg_pos=0.5596; %cg location (relative to the fuselage)
V_h=0.675; %horizontal tail volume coefficient
K_c=1.4; %tail arm correction factor

%% 2. Tail arm and surface computation.

%l_opt=K_c*sqrt((4*m_chord*S*V_h)/(pi*D_f)); (not use due to length limitation)
l_opt=10.2;
S_h=(m_chord*S*V_h)/l_opt;

%% 3. Lift needed for cruise flight stability

Cm_ow=Cm_o_w*(AR_w*cos(sweep_angle)^2)/(AR_w+2*cos(sweep_angle))+0.01*alpha_t;
X_apex=-h_o*m_chord+cg_pos*L_f+dist_Cg_Ca_w;
x_cg=h_o*m_chord-dist_Cg_Ca_w;
h=x_cg/m_chord;
Cl_h=(Cm_ow+Cl_w*(h-h_o))/V_h;

%% 4. Sizing calculations

AR_h=0.5*AR_w;
TR_h=TR_w;
alpha_twist_h=twist_angle;
sweep_angle_h=sweep_angle;
b_h=sqrt(AR_h*S_h);
C_h_mean=S_h/b_h;
C_h_root=C_h_mean*(3/2)*((1+TR_h)/(1+TR_h+TR_h^2));
C_h_tip=C_h_root*TR_h;

%% 5. Lift computation.

Cl_alpha_h=(0.6243*2)/(deg2rad(5)*2);
CL_alpha_h=Cl_alpha_h/(1+(Cl_alpha_h)/(pi*AR_h));
alpha_h=Cl_h/CL_alpha_h;
CL_tail=lift_calculator(S_h,AR_h,TR_h,alpha_h,CL_alpha_h,b_h,C_h_mean,C_h_root,alpha_twist_h);

%% 6. Iterative process

while abs(CL_tail-Cl_h)>0.002
if CL_tail-Cl_h<0
alpha_h=alpha_h+0.001;
elseif CL_tail-Cl_h>0
alpha_h=alpha_h-0.001;
end
CL_tail=lift_calculator(S_h,AR_h,TR_h,alpha_h,CL_alpha_h,b_h,C_h_mean,C_h_root,alpha_twist_h);
end

%% 7. Tail setting calculations.

alpha_w=deg2rad(1.75);
Cl_alpha_w=(0.6765+0.4404)/(deg2rad(5));
eps_o=(2*Cl_w)/(pi*AR_w);
der_eps=(2*Cl_alpha_w)/(pi*AR_w);
eps=eps_o+der_eps*alpha_w;
i_h=alpha_h+eps;

%% 8. Stability criteria

eta_h=l_opt*S_h/(S*m_chord);
C_m_alpha_h=Cl_alpha_w*(h-h_o)-CL_alpha_h*eta_h*(S_h/S);



