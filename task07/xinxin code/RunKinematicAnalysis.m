
clear all
clc
close all

%% Initial values of the system
l1 = 0.2;
l2 = 0.4;
l3 = 0.3;
d1 = 0.35;
d2 = 0.1;
phinia1 = 2.36;
phinia2 = 0.57;
phinia3 = 2.11;
dphinia1 = 6.28;

q0=[phinia1 phinia2 phinia3];
dq0=[dphinia1 0 0];
ddq0=[0 0 0];

%% Solution values
%maximum number of iterations
maxiter=100;
%parameter to establish the convergenvce criteria
epsilon=0.0000001;

cnt = 0;
tini = 0;  %Initial time
dt = 0.00001;  %Time step
tend = 1; %End time
% initial evaluation of constrains
t=0;
C = ConsEqua4Bar(q0,t);
Cq = Jaco4Bar(q0);

%% initial values of the system
q(:,1)=q0;
dq(:,1)=dq0;
ddq(:,1)=ddq0;
tsim(1)=tini; 

for t = tini:dt:tend
  
cnt = cnt+1; %Column counter to form the result matrix.
nloop=0;
nconverg=0;

while nloop<maxiter

%% Newton Difference
deltaq = -Cq^(-1)*C;
% Coordinate update
q(:,cnt)=q(:,cnt)+deltaq;

        nloop=nloop+1;
        nconverg=nconverg+1;    
 %%%% Converging criteria (error < epsilon)
maxdx=max(abs(deltaq));
        if maxdx<epsilon
            nloop=maxiter;
        end
 % constrains update       
tmp=q(:,cnt);
C = ConsEqua4Bar(tmp,t);
% % jacobian update
Cq = Jaco4Bar(tmp);
tsim(cnt+1)=t;
end

%% Estimation of positions for next step
q(:,cnt+1)=q(:,cnt);
%Presentation and store of the result values
q;
nconverg;
%Determination of the velocities (dq)
Ct = ConsTime(tmp,t);
% Velocity  
dq(:,cnt) = Cq^(-1)*-Ct;
dq(:,cnt+1) = dq(:,cnt);

%% Determination of the accelerations
tmpdq = dq(:,cnt);
Ctt = ConsTimeTime(tmpdq,t); 
Cqqq = Consqq(tmp,tmpdq,t);
Cqt = Consqt(tmpdq,t);
 
ddq(:,cnt)=Cq^(-1)*(-Ctt-Cqqq-Cqt);
ddq(:,cnt+1)=dq(:,cnt);
end

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 8 6]);
plot(tsim,dq(1,:),'LineWidth', 1.5)
hold on
plot(tsim,dq(2,:),'LineWidth', 1.5)
hold on
plot(tsim,dq(3,:),'LineWidth', 1.5)
legend({'input link','coupler','output link'},'FontUnits','points','interpreter','latex','FontSize',8,'FontName','Times','Location','NorthWest')
xlabel('Time s','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
ylabel('Angular Velocity rad/s','FontUnits','points','interpreter','latex','FontSize',11,'FontName','Times')
print -depsc2 angularvelocity.eps
grid on
hold off



% % Plottings and animation
%  
% q = q';   % Positions
% dq = dq';  % Velocities.
% ddq = ddq'; % Accelerations. 
% 
% % Positions of the end points 
% pa = [zeros(1,length(q)) zeros(1,length(q))];
% pb = [l1*cos(q(:,1)) l1*sin(q(:,1))];
% pb2 = [l1*cos(q(:,1))+l2*cos(q(:,2)) l1*sin(q(:,1))+l2*sin(q(:,2))];
% pc = [d1*ones(1,length(q)) d2*ones(1,length(q))];
% 
% hold off
% 
% % Visualization of numerical results
% disp('Press any key for visualization')
% pause
% 
% for i=1:max(size(pa))
% 	plot([pa(i,1) pb(i,1) pb2(i,1) pc(i,1)],[pa(i,2) pb(i,2) pb2(i,2) pc(i,2)],'b-',-.1,-1,'wo',2.5,1,'wo',0,0,'ro')
% 	pause(.005)
% end
% 
% hold off