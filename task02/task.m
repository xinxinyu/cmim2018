%clear all
format long
p.m=1;
p.k=4;
p.c=0.1;
p.A=100;
p.omga=20;
tic
v0=[-1 0]; %% Initial displacement and velocity of the system
dt = 0.02;
TT=20;
% options = odeset('MaxStep',1e-5);
% %options=[];
% [t,v]=ode45(@(t,v) SysPosition(t,v,p),[0 : dt : TT], v0);
% % [t,v]=ode15s(@SysPosition,(0 : dt : TT), v0,options, p);
% % [t,v]=ode15s(@(t,v) SysPosition(t,v,p),[0 : dt : TT], v0);
% % [t,v]=ode23(@(t,v) SysPosition(t,v,p),[0 : dt : TT], v0);
% v
RelTol=1e-6;
AbsTol=1e-9; 
options = odeset('AbsTol',AbsTol,'RelTol',RelTol);
[t1,v1]=ode45(@(t,v) SysPosition(t,v,p),[0 : dt : TT], v0,options);
% [t2,v2]=ode15s(@(t,v) SysPosition(t,v,p),[0 : dt : TT], v0,options);
% [t3,v3]=ode23(@(t,v) SysPosition(t,v,p),[0 : dt : TT], v0,options);

plot(v1(:,1),'LineWidth', 1.5 );

hold on
plot(v1(:,2),'--','LineWidth', 1.5 );
hold on
legend('position','velocity')
% plot(v1(:,1),'LineWidth', 1.5 );
% hold on
% plot(v2(:,2),'--','LineWidth', 1.5 );
% hold on
% plot(v2(:,1),'--','LineWidth', 1.5 );
% hold on
% plot(v3(:,2),':','LineWidth', 1.5 );
% plot(v3(:,1),'--','LineWidth', 1.5 );
% hold on
% hold off
grid on
% 
% legend('velocity-ode45','position-ode45','velocity-ode15s','position-ode15s','velocity-ode23','position-ode23')
toc

% figure(1)
% t=0:dt:TT;
% plot(t,v(:,1))
% hold on 
% plot(t,v(:,2))
% xlabel('Time [s]')
% legend('positon','velocity')
% hold off

