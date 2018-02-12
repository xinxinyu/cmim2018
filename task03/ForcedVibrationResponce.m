%% Force vibration responce
%written by Xinxin
% 08.02.2018


%% Define symbolic variables
clc
clear
close all

syms  K F0 X0 A OMEGA omega zeta  f real
% assumeAlso(m > 0);
% assumeAlso(k > 0);
% assumeAlso(c > 0);

syms x(t)
%% Analyze the system with the initial condition

xp = diff(x, t);
xpp = diff(xp, t);

% omega=sqrt(k/m);zeta=c/2/sqrt(k*m); K=1/k
Eq1 = (1/omega)^2*xpp+(2*zeta/omega)*xp+x-K*F0*sin(OMEGA*t);

solution_damping = dsolve(Eq1 == 0, [xp(0) == 0, x(0) == 0]);
symvar(solution_damping)
disp(solution_damping);


solution_damping_rewritten=simplify(expand(solution_damping));
disp(solution_damping_rewritten);
symvar(solution_damping_rewritten);

fplot(subs(solution_damping_rewritten, [ F0, K, OMEGA, omega, zeta], [10,1,6,2,0.5]), [0, 10])

%% Analyze Force vibration responce
Eq2=K*F0/sqrt((1-(omega/OMEGA)^2).^2+(2*zeta*(omega/OMEGA)).^2) -X0;
solution_amplititude= solve(Eq2==0,X0);
disp(solution_amplititude);
symvar(solution_amplititude);

solution_amplititude_rewrite = subs(solution_amplititude, omega/OMEGA, f);
disp(solution_amplititude_rewrite) 
symvar(solution_amplititude_rewrite) 

solution_amplititude_rewritten = simplify(solution_amplititude_rewrite/(F0*K));
disp(solution_amplititude_rewritten) 
symvar(solution_amplititude_rewritten) 

%% First method to Plot the figure

figure
fplot(subs(solution_amplititude_rewritten, [zeta], [0:0.1:1]), [0, 10],'LineWidth', 1.5)
grid on

xlim([0,2.5]);
ylim([0,5]);
xlabel('Frequency Ratio $$r={f\over f_{n}}$$','FontUnits','points','interpreter','latex','FontSize',12,'FontName','Times')
ylabel('Amplication Ratio $X\frac{k}{F_{0}}$  ','FontUnits','points','interpreter','latex','FontSize',12,'FontName','Times')
title('Amplitude','FontUnits','points','FontWeight','normal','FontSize',16,'FontName','Times')
grid on
print -depsc2 myplot.eps
print -dpdf myplot.pdf
hold off


%% Second method to Plot the figure
i=1;
for zeta=0:0.125:1.125
    i=i+1;
M(i)=K*F0/sqrt((1-(omega*t/OMEGA)^2).^2+(2*zeta(1,:)*(omega*t/OMEGA)).^2);
% M(i)=f/sqrt(((1-(OMEGA*t/omega)^2).^2+(2*zeta(1,:)*(OMEGA*t/omega).^2)));
end 
disp(M) 
symvar(M) 
fplot(subs(M, [ F0, K, OMEGA, omega], [1,1,6,2]), [0, 5],'LineWidth', 1.0)
grid on






%% Second method to plot the figure with loop
% % wn=2;   %%wn=sqrt(k/m)
% % wdr=3;
% fo=1;
% % f=linspace(0,4,100);
% zeta1=0;
% n=1;
% for i=0:0.025:2.5
%     A1(n)=fo/sqrt((1-i^2).^2+(2*zeta1*i).^2);
% 
%     n=n+1;
%     disp(i)
% end
% 
% zeta2=0.1;
% n=1;
% for i=0:0.025:2.5
%     A2(n)=fo/sqrt((1-i^2).^2+(2*zeta2*i).^2);
% 
%     n=n+1;
% end
% 
% zeta3=0.2;
% n=1;
% for i=0:0.025:2.5
%     A3(n)=fo/sqrt((1-i^2).^2+(2*zeta3*i).^2);
% 
%     n=n+1;
% end
% 
% zeta4=0.3;
% n=1;
% for i=0:0.025:2.5
%     A4(n)=fo/sqrt((1-i^2).^2+(2*zeta4*i).^2);
% 
%     n=n+1;
% end
% 
% zeta5=0.4;
% n=1;
% for i=0:0.025:2.5
%     A5(n)=fo/sqrt((1-i^2).^2+(2*zeta5*i).^2);
% 
%     n=n+1;
% end
% 
% zeta6=0.5;
% n=1;
% for i=0:0.025:2.5
%     A6(n)=fo/sqrt((1-i^2).^2+(2*zeta6*i).^2);
% 
%     n=n+1;
% end
% 
% zeta7=1;
% n=1;
% for i=0:0.025:2.5
%     A7(n)=fo/sqrt((1-i^2).^2+(2*zeta7*i).^2);
%     n=n+1;
% end
% figure('Units','inches','Position',[0 0 6 4],'PaperPositionMode','auto');
% plot(0:0.025:2.5,A1,'LineWidth', 1.5)
% hold on
% plot(0:0.025:2.5,A2,'LineWidth', 1.5)
% hold on
% plot(0:0.025:2.5,A3,'LineWidth', 1.5)
% hold on
% plot(0:0.025:2.5,A4,'LineWidth', 1.5)
% hold on
% plot(0:0.025:2.5,A5,'LineWidth', 1.5)
% hold on
% plot(0:0.025:2.5,A6,'LineWidth', 1.5)
% hold on
% plot(0:0.025:2.5,A7,'LineWidth', 1.5)
% xlim([0,2.5]);
% ylim([0,5]);
% hold off

% xlabel('Frequency Ratio $$r={f\over f_{n}}$$','FontUnits','points','interpreter','latex','FontSize',12,'FontName','Times')
% ylabel('Amplication Ratio $X\frac{k}{F_{0}}$  ','FontUnits','points','interpreter','latex','FontSize',12,'FontName','Times')
% % legend({'2'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Times','Location','NorthEast')
% title('Amplitude','FontUnits','points','FontWeight','normal','FontSize',16,'FontName','Times')
% grid on
% print -depsc2 myplot.eps
% print -dpdf myplot.pdf
% hold off

