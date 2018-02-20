%% test mass spring system

m=1;
k=100;

dt = 1e-6;
tk = 5;

y0=[1 0];

[T, Y] = odeFE(@(t,y) Spring(t,y,m,k), [dt, tk], y0);


V=0.5*k*Y(1,:).^2;
U=0.5*m*Y(2,:).^2;
Total=V+U;

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8 6]);
plot(0:dt:tk,Y(1,:),'LineWidth', 1.2)
hold on
plot(0:dt:tk,Total,'LineWidth', 1.2)
xlabel('Time','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
ylabel('Diaplacement-Energy','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
legend({'position','total energy'},'FontUnits','points','interpreter','latex','FontSize',6,'FontName','Times','Location','NorthWest')
grid on
print -depsc2 task3.eps
hold off