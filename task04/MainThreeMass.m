clear all
clc;

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
%% Calculate egien frequency
e = eig(K,M);
ee =sqrt(e);

%% initial condition of the system
y0=[1 1.5 2 1.5 3 4];
%% Calculate the displacement and velocity of the system

dt = 10e-5;
tk = 0.5;
RelTol=1e-6;
AbsTol=1e-9; 
options = odeset('AbsTol',AbsTol,'RelTol',RelTol,'Stats','on');
%%%%%%%%%%%%%%%%%%%%%%%%%%%% ode45 ****************
tic
[T, Y1] = ode45(@(t,y) ThreeSpring(t,y,M,K), [0:dt:tk], y0,options);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%% ode15s ****************
tic
[T, Y2] = ode15s(@(t,y) ThreeSpring(t,y,M,K), [0:dt:tk], y0,options);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%% odeRK4 ****************
tic
[T, Y3] = odeRK4(@(t,y) ThreeSpring(t,y,M,K), [dt tk], y0);
toc
Y3=Y3';
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Semi- implicit Euler Method ****************
fun_f = @(t,v) v;
fun_g = @(t,u) -1*(M^-1)*K*u;
u0 = y0(1:3); % initial displacement
v0 = y0(4:6); % initial velocity
% [T,u,v] = simEuler(fun_f,fun_g,[0:dt:tk],u0,v0);
tic
[T,u,v] = odeSemiFE(fun_f,fun_g,[0:dt:tk],u0,v0);
toc
Y4 = [u;v];
Y4 = Y4';
%%%%%%%%%%%%%%%%%%%%%%%%%%%% BackwardEuler Method ****************
[T,Y] = BackwardEuler(@(t,y) ThreeSpring(t,y,M,K),[0:dt:tk],y0);

%% Calculate the Modal Coordinates

[T, Y_tilda1] = ode45(@(t,y) modal_coordinates(t,y,M,K), [0:dt:tk], y0);
[T, Y_tilda2] = ode15s(@(t,y) modal_coordinates(t,y,M,K), [0:dt:tk], y0);
[T, Y_tilda3] = odeRK4(@(t,y) modal_coordinates(t,y,M,K), [dt tk], y0);
Y_tilda3 = Y_tilda3';


%% Calculate accuracy with four methods
rmse_ode15 = sqrt(sum(Y2(:,1)-Y1(:,1)).^2/length(Y1(:,1)));
rmse_odeRK4 = sqrt(sum(Y3(:,1)-Y1(:,1)).^2/length(Y1(:,1)));
rmse_odeSemiFE = sqrt(sum(Y4(:,1)-Y1(:,1)).^2/length(Y1(:,1)));
gap1=Y2(:,1)-Y1(:,1);
gap2=Y3(:,1)-Y1(:,1);
gap3=Y4(:,1)-Y1(:,1);
errorV1=Y2(:,4)-Y1(:,4);
errorV2=Y3(:,4)-Y1(:,4);
errorV3=Y4(:,4)-Y1(:,4);

%% Plot all the figures
%%Plot the accuracy of displacement the system

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8 6]);
plot(0:dt:tk,gap1,'LineWidth',1.2,'color','black')
hold on
plot(0:dt:tk,gap2,'--','LineWidth', 1.2,'color','black')
hold on
plot(0:dt:tk,gap3,':','LineWidth', 1.2,'color','black')
xlabel('Time (s)','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
ylabel('error (m)','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
legend({'ode45-ode15s','ode45-odeRK4','ode45-odeSIE'},'FontUnits','points','interpreter','latex','FontSize',6,'FontName','Times','Location','NorthWest')
grid on
print -depsc2 errorofdisplacement.eps
hold off

%%Plot the accuracy of velocity the system

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8 6]);
plot(0:dt:tk,errorV1,'LineWidth', 1.2,'color','black')
hold on
plot(0:dt:tk,errorV2,'--','LineWidth', 1.2,'color','black')
hold on
plot(0:dt:tk,errorV3,':','LineWidth', 1.2,'color','black')
xlabel('Time (s)','FontUnits','points','interpreter','latex','FontSize',10,'FontName','Times')
ylabel('error (m/s)','FontUnits','points','interpreter','latex','FontSize',10,'FontName','Times')
legend({'ode45-ode15s','ode45-odeRK4','ode45-odeSIE'},'FontUnits','points','interpreter','latex','FontSize',6,'FontName','Times','Location','NorthWest')
grid on
print -depsc2 errorofvelocity.eps
hold off


%%Plot the difference of displacement between general coordinate with modal coordinates

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8 6]);
plot(0:dt:tk,Y1(:,1),'LineWidth', 1.2, 'color','black')
xlabel('Time (s)','FontUnits','points','interpreter','latex','FontSize',10,'FontName','Times')
ylabel('General displacement (m)','FontUnits','points','interpreter','latex','FontSize',10,'FontName','Times')
grid on
print -depsc2 Generaldisplacement.eps
hold off

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8 6]);
plot(0:dt:tk,Y_tilda1(:,1),'LineWidth', 1.2,'color','black')
xlabel('Time (s)','FontUnits','points','interpreter','latex','FontSize',10,'FontName','Times')
ylabel('Modal displacement (m)','FontUnits','points','interpreter','latex','FontSize',10,'FontName','Times')
grid on
print -depsc2 Modaldisplacement.eps
hold off


set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8 6]);
plot(0:dt:tk,Y1(:,4),'LineWidth', 1.2, 'color','black')
xlabel('Time (s)','FontUnits','points','interpreter','latex','FontSize',10,'FontName','Times')
ylabel('General velocity (m/s)','FontUnits','points','interpreter','latex','FontSize',10,'FontName','Times')
grid on
print -depsc2 Generalvelocity.eps
hold off

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8 6]);
plot(0:dt:tk,Y_tilda1(:,4),'LineWidth', 1.2,'color','black')
xlabel('Time (s)','FontUnits','points','interpreter','latex','FontSize',10,'FontName','Times')
ylabel('Modal velocity (m/s)','FontUnits','points','interpreter','latex','FontSize',10,'FontName','Times')
grid on
print -depsc2 Modalvelocity.eps
hold off






set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8 6]);
plot(0:dt:tk,Y1(:,1),'LineWidth', 1.2)
hold on
plot(0:dt:tk,Y_tilda1(:,1),'LineWidth', 1.2)
hold on
xlabel('Time','FontUnits','points','interpreter','latex','FontSize',10,'FontName','Times')
ylabel('error','FontUnits','points','interpreter','latex','FontSize',10,'FontName','Times')
legend({'General displacement','Modal displacement'},'FontUnits','points','interpreter','latex','FontSize',6,'FontName','Times','Location','NorthWest')
grid on
print -depsc2 tildalofDisplacement.eps
hold off

%%Plot the difference of velocity between general coordinate with modal coordinates

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8 6]);
plot(0:dt:tk,Y1(:,4),'LineWidth', 1.2)
hold on
plot(0:dt:tk,Y_tilda1(:,4),'LineWidth', 1.2)
hold on
xlabel('Time','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
ylabel('error','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
legend({'General velocity','Modal velocity'},'FontUnits','points','interpreter','latex','FontSize',6,'FontName','Times','Location','NorthWest')
grid on
print -depsc2 tildalofVelocity.eps
hold off

% 

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8 6]);
plot(0:dt:tk,Y1(:,1),'LineWidth', 1.5)
hold on
plot(0:dt:tk,Y1(:,2),'LineWidth', 1.5)
hold on
plot(0:dt:tk,Y1(:,3),'LineWidth', 1.5)
% hold on
% plot(Y4(:,1))

xlabel('Time','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
ylabel('Position','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
legend({'Mass 1','Mass 2','Mass 3'},'FontUnits','points','interpreter','latex','FontSize',8,'FontName','Times','Location','NorthEast')
grid on
print -depsc2 myplot.eps
print -dpdf myplot.pdf
hold off

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8 6]);
plot(0:dt:tk,Y1(:,4),'LineWidth', 1.2)
hold on
plot(0:dt:tk,Y1(:,5),'LineWidth', 1.2)
hold on
plot(0:dt:tk,Y1(:,6),'LineWidth', 1.2)
% hold on
% plot(Y4(:,1))
xlabel('Time','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
ylabel('Velocity','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
legend({'Mass 1','Mass 2','Mass 3'},'FontUnits','points','interpreter','latex','FontSize',9,'FontName','Times','Location','NorthEast')
grid on
print -depsc2 myplot1.eps
print -dpdf myplot1.pdf
hold off
