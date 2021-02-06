
clc;
clear all;
close all;

C_L = [-2:0.01:3];

C_D_clean = 0.020+0.0438*C_L.^2;
C_D_to_flaps = 0.035+0.0448*C_L.^2;
C_D_l_flaps = 0.085+0.0517*C_L.^2;
C_D_l_gear = 0.040+0.0438*C_L.^2;

fig1=figure(1);
set(fig1,'Renderer', 'painters', 'Position', [400 400 500 350]);
hold on

grid on;
grid minor;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.2;
plot(C_L,C_D_clean,'DisplayName','$C_D$ (Clean)');
plot(C_L,C_D_to_flaps,'DisplayName','$C_D$ (Take-off flaps)');
plot(C_L,C_D_l_flaps,'DisplayName','$C_D$ (Landing flaps)');
plot(C_L,C_D_l_gear,'DisplayName','$C_D$ (Gear)');
hold off;

% Axis format
set(gca,'TickLabelInterpreter','latex','fontsize',12)
xlabel('Lift coefficient $C_L$','interpreter','latex','FontSize',14)
ylabel('Drag coefficient $C_D$','interpreter','latex','FontSize',14)
legend('location','northwest','interpreter','latex');

print(gcf,'drag_polar.png','-dpng','-r800');

