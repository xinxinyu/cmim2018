
% clear all
% clc;
dt = 10e-6;
tk = 1;

m1=2;
m2=1.5;
m3=1;
M=diag([m1,m2,m3]);
k1=20000;
k2=15000;
k3=10000;
K=[k1+k2 -k2 0;
    -k2 k2+k3 -k3;
    0 -k3 k3];

[psi,lambda]=eig(K,M);

y0=[0 0 3 0 0 0];

% M_tilda1 = psi(:,1)'*M*psi(:,1);
% K_tilda1 = psi(:,1)'*K*psi(:,1);
% M_tilda2 = psi(:,2)'*M*psi(:,2);
% K_tilda2 = psi(:,2)'*K*psi(:,2);
% M_tilda3 = psi(:,3)'*M*psi(:,3);
% K_tilda3 = psi(:,3)'*K*psi(:,3);
% 
% 
% [psi_t,lambda_t]=eig(K_tilda,M_tilda);

RelTol=1e-6;
AbsTol=1e-9; 
options = odeset('AbsTol',AbsTol,'RelTol',RelTol);

[T, Y1] = ode45(@(t,y) ThreeSpring(t,y,M,K), [0:dt:tk], y0,options);

[T, Y2] = ode15s(@(t,y) ThreeSpring(t,y,M,K), [0:dt:tk], y0,options);

[T, Y3] = odeRK4(@(t,y) ThreeSpring(t,y,M,K), [dt tk], y0);
Y3=Y3';

%%Semi- implicit Euler Method
fun_f = @(t,v) v;
fun_g = @(t,u) -1*(M^-1)*K*u;
u0 = y0(1:3);
v0 = y0(4:6);
% [T,u,v] = simEuler(fun_f,fun_g,[0:dt:tk],u0,v0);
% % % [T,u,v] = odeSemiFE(fun_f,fun_g,[0:dt:tk],u0,v0);
Y4 = [u;v];
Y4=Y4';
% [T, Y4] = odeSemiFE(@(t,y) ThreeSpring(t,y,M,K), [dt tk], y0,K,M);
% Y4=Y4';

rmse = @(x1, x2) sqrt(sum((x1-x2).^2) / length(x1));

disp(['Time step: ',num2str(dt), ' norm(err) = ', num2str(rmse(Y1(:,1),Y2(:,1)))])
disp(['Time step: ',num2str(dt), ' norm(err) = ', num2str(rmse(Y1(:,1),Y3(:,3)))])
disp(['Time step: ',num2str(dt), ' norm(err) = ', num2str(rmse(Y2(:,1),Y3(:,3)))])


% gap1=Y1(:,1)-Y2(:,1);
% gap2=Y1(:,1)-Y3(:,1);
% gap3=Y2(:,1)-Y3(:,1);
% 
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 8 6]);
% plot(0:dt:tk,gap1,'LineWidth', 1.5)
% hold on
% plot(0:dt:tk,gap2,'LineWidth', 1.5)
% hold on
% plot(0:dt:tk,gap3,'LineWidth', 1.5)
% % hold on
% % plot(Y4(:,1))
% 
% xlabel('Time','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
% ylabel('error','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
% legend({'ode45-ode15s','ode45-odeRK4','ode15s-odeRK4'},'FontUnits','points','interpreter','latex','FontSize',8,'FontName','Times','Location','NorthWest')
% grid on
% print -depsc2 myplot3.eps
% print -dpdf myplot3.pdf
% hold off
% 
% 
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 8 6]);
% plot(0:dt:tk,Y1(:,1),'LineWidth', 1.5)
% hold on
% plot(0:dt:tk,Y1(:,2),'LineWidth', 1.5)
% hold on
% plot(0:dt:tk,Y1(:,3),'LineWidth', 1.5)
% % hold on
% % plot(Y4(:,1))
% 
% xlabel('Time','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
% ylabel('Position','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
% legend({'Mass 1','Mass 2','Mass 3'},'FontUnits','points','interpreter','latex','FontSize',8,'FontName','Times','Location','NorthEast')
% grid on
% print -depsc2 myplot.eps
% print -dpdf myplot.pdf
% hold off
% 
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 8 6]);
% plot(0:dt:tk,Y1(:,4),'LineWidth', 1.2)
% hold on
% plot(0:dt:tk,Y1(:,5),'LineWidth', 1.2)
% hold on
% plot(0:dt:tk,Y1(:,6),'LineWidth', 1.2)
% % hold on
% % plot(Y4(:,1))
% 
% xlabel('Time','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
% ylabel('Velocity','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
% legend({'Mass 1','Mass 2','Mass 3'},'FontUnits','points','interpreter','latex','FontSize',9,'FontName','Times','Location','NorthEast')
% grid on
% print -depsc2 myplot1.eps
% print -dpdf myplot1.pdf
% hold off
