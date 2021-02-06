function man_diagr_density(rho_SL,MTOW,Mach_cruise,S_w,C_Nmax,C_NmaxTO,V_cruise,T,R,gamma,i)
% Function that graphs the manoeuvring diagram and determines its
% velocities according to the density selected

%% Maximum load factor (n) conditions
x_velocity_vector = [20:1:500];
n_max_lim = 3.8;
n_min_lim = -2;
n_max = 2.1 + 24000/(MTOW/0.453592 + 10000);
if n_max < 2.5
    n_max = 2.5;
end

%% Diving speed (velocidad de proyecto de picado)
M_D = 0.07 + Mach_cruise;
V_D = sqrt(gamma*R*T)*M_D; % [m/s]

%% Stall curve C_Nmax
V_s_C_Nmax=@(n) sqrt(n*MTOW/(0.5*rho_SL*S_w*C_Nmax));
V_s1=double(subs(V_s_C_Nmax,1));
V_A=double(subs(V_s_C_Nmax,n_max));

%% Stall curve C_NmaxTO
n_max_F=2; % With flaps
V_s_C_NmaxTO=@(n) sqrt(n*MTOW/(0.5*rho_SL*S_w*C_NmaxTO));
V_s1_TO=double(subs(V_s_C_NmaxTO,1));
V_s_F=double(subs(V_s_C_NmaxTO,n_max_F));
V_A_F=double(subs(V_s_C_NmaxTO,n_max_F));
V_F=1.6*V_s1_TO;
 
%% Stall curve n_min
n_min=-1;
V_s_n_min=@(n) sqrt(n*MTOW/(0.5*rho_SL*S_w*(-0.8)*C_Nmax));
V_s_minus_1=double(subs(V_s_n_min,n_min));
if i==2
     figure(i)
elseif i==3
     openfig('SL_man_diag')
end

hold on;
vec_n_V_stall=linspace(0,n_max,100);
plot(double(subs(V_s_C_Nmax,vec_n_V_stall)),vec_n_V_stall,'b');

vec_n_V_stall_F=linspace(0,n_max_F,100);
plot(double(subs(V_s_C_NmaxTO,vec_n_V_stall_F)),vec_n_V_stall_F,'b');

vec_V_n_max_F=linspace(V_A_F,V_F,2);
plot(vec_V_n_max_F,[n_max_F,n_max_F],'b');

vec_V_n_max=linspace(V_A,V_D,2);
plot(vec_V_n_max,[n_max,n_max],'b');
 
vec_n_V_D=linspace(0,n_max,2);
plot([V_D, V_D],vec_n_V_D,'b');

vec_n_neg_V_stall=linspace(0,n_min,100);
plot(double(subs(V_s_n_min,vec_n_neg_V_stall)),vec_n_neg_V_stall,'b');

vec_V_n_min=linspace(V_s_minus_1,V_cruise,2);
plot(vec_V_n_min,[n_min,n_min],'b');

plot([V_cruise V_D], [n_min 0],'b');

if i==2 % Sea Level case:
    savefig('SL_man_diag.fig');
end

ylabel("Load factor $n$");
xlabel("Velocity $\left[ \mathrm{m}/\mathrm{s} \right]$");

if i==2 % Sea Level case:
    title("Manoeuvering diagram (SL)");
elseif i==3 % Cruise case:
    title("Manoeuvering diagram (h=12km)");
end

xlim([0 300]);
ylim([-1.5 3]);
xticks([0 100 200 300])
xticklabels([0,100,200,300])
grid on;
grid minor;
box on;

end