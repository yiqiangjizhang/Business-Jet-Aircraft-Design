function Cl_h=lift_calculator(S_h,AR_h,TR_h,alpha_h,CL_alpha_h,b_h,C_h_mean,C_h_root,alpha_twist_h)
N = 9; % (number of segments-1)
alpha_twist=rad2deg(alpha_twist_h); % Twist angle (deg)
a_h=rad2deg(alpha_h); % tail angle of attack (deg)
alpha_0=rad2deg(0); % zero-lift angle of attack (deg)
theta = pi/(2*N):pi/(2*N):pi/2;
alpha=a_h+alpha_twist:-alpha_twist/(N-1):a_h;
% segments angle of attack
z = (b_h/2)*cos(theta);
c = C_h_root * (1 - (1-TR_h)*cos(theta)); % Mean
%Aerodynamics chord at each segment
mu = c * CL_alpha_h / (4 * b_h);
LHS = mu .* (alpha-alpha_0)/57.3; % Left Hand Side
% Solving N equations to find coefficients A(i):
for i=1:N
for j=1:N
B(i,j)=sin((2*j-1) * theta(i)) * (1+(mu(i) *(2*j-1))/sin(theta(i)));
end
end
A=B\transpose(LHS);
for i = 1:N
sum1(i) = 0;
sum2(i) = 0;
for j = 1 : N
sum1(i) = sum1(i) + (2*j-1) * A(j)*sin((2*j-1)*theta(i));
sum2(i) = sum2(i) + A(j)*sin((2*j-1)*theta(i));
end
end
Cl_h= pi * AR_h * A(1);