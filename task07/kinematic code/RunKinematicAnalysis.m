%% written by xinxin 21/3/2018
clear all
clc
close all

%% initial data
% model definition

syms l1 l2 l3 l4 real
syms x1 y1 phi1 x2 y2 phi2 x3 y3 phi3 x4 y4 phi4 d1 d2 real


d1 = 2;
d2 = 0;
body1.type = 'slender_link';
body1.l = 2;
body1.uix = -1;
body1.uiy = 0;
body1.ujx = 1;
body1.ujy = 0;
body1.driver = 1;
body1.driver_position = pi/2;
body1.driver_speed = 3.14;

body2.type = 'slender_link';
body2.l = 2;
body2.uix = -1;
body2.uiy = 0;
body2.ujx = 1;
body2.ujy = 0;
body2.driver = 0;
body2.driver_position = 0;
body2.driver_speed = 0;

body3.type = 'slender_link';
body3.l = 2;
body3.uix = -1;
body3.uiy = 0;
body3.ujx = 1;
body3.ujy = 0;
body3.driver = 0;
body3.driver_position = 0;
body3.driver_speed = 0;

model_def.bodies = [body1, body2, body3];

joint1.type = 'revolute';
joint1.ibody = 0;
joint1.jbody = 1;

joint2.type = 'revolute';
joint2.ibody = 1;
joint2.jbody = 2;

joint3.type = 'revolute';
joint3.ibody = 2;
joint3.jbody = 3;

joint4.type = 'revolute';
joint4.ibody = 3;
joint4.jbody = 0;

model_def.joints = [joint1, joint2, joint3, joint4];


q0 = [0 1 pi/2 1 2 0 2 1 -(pi/2)];
dq0 = [3.14 0 0 0 0 0 0 0 0];
ddq0 = zeros (1,length(q0));

tini = 0;  %Initial time
dt = 0.001;  %Time step
tend = 0.5; %End time
[t,q,dq,ddq] = kinematic_solution(model_def,tini:dt:tend,q0,dq0,ddq0,d1,d2);


%% input the data into Matlab
[A,b]=xlsread('Z:\cmim2018\cmim2018\task07\kinematic code\input.xlsx');
c=str2double(b);
d=c';   

%% 
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 12 10]);
xlim([0.176 0.192]);
ylim([0.82 0.86])
plot(tini:dt:tend,q(2,2:end),'LineWidth', 1.2,'color','black')
hold on
plot(tini:dt:tend,d(2,1:end),':','LineWidth', 1.2,'color','black')
% % hold on
% % plot(tini:dt:tend,q(8,2:end),'LineWidth', 1.5)
legend({'Matlab','Adams','output link'},'FontUnits','points','interpreter','latex','FontSize',8,'FontName','Times')
xlabel('Time (s)','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
ylabel('position of the driving bar (m)','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
grid on
print -depsc2 comparecut.eps
hold off


