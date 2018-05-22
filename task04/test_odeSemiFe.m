%% test mass spring system
clear up;
clc;
m=1;
k=100;

dt = 1e-6;
tk = 2.5;

y0=[1 0];
fun_f = @(t,v) v;
fun_g = @(t,u) -1*(m^-1)*k*u;
u0 = y0(1); % initial displacement
v0 = y0(2); % initial velocity
% [T,u,v] = simEuler(fun_f,fun_g,[0:dt:tk],u0,v0);
[T,u,v] = odeSemiFE(fun_f,fun_g,[0:dt:tk],u0,v0);
Y4 = [u;v];

[T, Y1] = odeFE(@(t,y) Spring(t,y,m,k), [dt, tk], y0);
[T, Y2] = odeRK4(@(t,y) Spring(t,y,m,k), [dt tk], y0);

gap1=Y2(1,:)-Y1(1,:);
gap2=Y2(1,:)-Y4(1,:);

V=0.5*k*Y(1,:).^2;
U=0.5*m*Y(2,:).^2;
Total=V+U;

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8 6]);
plot(0:dt:tk,Y4(1,:),'LineWidth',1.2,'color','black')
hold on
plot(0:dt:tk,Y4(2,:),'--','LineWidth', 1.2,'color','black')
xlabel('Time (s)','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
ylabel('Position (m)-Velocity (m/s)','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
legend({'position','velocity'},'FontUnits','points','interpreter','latex','FontSize',6,'FontName','Times','Location','NorthWest')
grid on
print -depsc2 task4.eps
hold off



set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8 6]);
plot(0:dt:tk,gap1(1,:),'LineWidth', 1.2,'color','black')
hold on
plot(0:dt:tk,gap2(1,:),'--','LineWidth', 1.2,'color','black')
hold on
xlabel('Time (s)','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
ylabel('Error between different solvers','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
legend({'Error between Runge-Kutta and Forward Euler','Error between Runge-Kutta and Semi-implicit Euler'},'FontUnits','points','interpreter','latex','FontSize',6,'FontName','Times','Location','NorthWest')
grid on
print -depsc2 task5_error.eps
hold off


