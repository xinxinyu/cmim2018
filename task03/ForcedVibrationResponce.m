clear all

%% Symbolic toolbox tutorial - in appliction to SDOF systems
% % Open help for symbolic toolbox.

%% Define symbolic variables
%%
clc
clear
close all

syms  K f A OMEGA omega zeta real
% assumeAlso(m > 0);
% assumeAlso(k > 0);
% assumeAlso(c > 0);

syms x(t)
%% Analyze damping vibrations

xp = diff(x, t);
xpp = diff(xp, t);

% omega=sqrt(k/m);zeta=c/2/sqrt(k*m); K=1/k
Eq = (1/omega)^2*xpp+(2*zeta/omega)*xp+x-K*f*sin(OMEGA*t);

solution_damping = dsolve(Eq == 0, [xp(0) == 0, x(0) == A]);
symvar(solution_damping)
disp(solution_damping);


solution_damping_rewritten=simplify(expand(solution_damping));
disp(solution_damping_rewritten);
symvar(solution_damping_rewritten);
% wn=2;   %%wn=sqrt(k/m)
% wdr=3;
fo=1;
% f=linspace(0,4,100);
zeta1=0;
n=1;
for i=0:0.025:2.5
    A1(n)=fo/sqrt((1-i^2).^2+(2*zeta1*i).^2);

    n=n+1;
    disp(i)
end

zeta2=0.1;
n=1;
for i=0:0.025:2.5
    A2(n)=fo/sqrt((1-i^2).^2+(2*zeta2*i).^2);

    n=n+1;
end

zeta3=0.2;
n=1;
for i=0:0.025:2.5
    A3(n)=fo/sqrt((1-i^2).^2+(2*zeta3*i).^2);

    n=n+1;
end

zeta4=0.3;
n=1;
for i=0:0.025:2.5
    A4(n)=fo/sqrt((1-i^2).^2+(2*zeta4*i).^2);

    n=n+1;
end

zeta5=0.4;
n=1;
for i=0:0.025:2.5
    A5(n)=fo/sqrt((1-i^2).^2+(2*zeta5*i).^2);

    n=n+1;
end

zeta6=0.5;
n=1;
for i=0:0.025:2.5
    A6(n)=fo/sqrt((1-i^2).^2+(2*zeta6*i).^2);

    n=n+1;
end

zeta7=1;
n=1;
for i=0:0.025:2.5
    A7(n)=fo/sqrt((1-i^2).^2+(2*zeta7*i).^2);
    n=n+1;
end
figure('Units','inches','Position',[0 0 6 4],'PaperPositionMode','auto');
plot(0:0.025:2.5,A1,'LineWidth', 1.5)
hold on
plot(0:0.025:2.5,A2,'LineWidth', 1.5)
hold on
plot(0:0.025:2.5,A3,'LineWidth', 1.5)
hold on
plot(0:0.025:2.5,A4,'LineWidth', 1.5)
hold on
plot(0:0.025:2.5,A5,'LineWidth', 1.5)
hold on
plot(0:0.025:2.5,A6,'LineWidth', 1.5)
hold on
plot(0:0.025:2.5,A7,'LineWidth', 1.5)
xlim([0,2.5]);
ylim([0,5]);
hold off
ax.XTickLabel = {'0','\pi','2\pi','3\pi'};
xlabel('Frequency Ratio $$r={f\over f_{n}}$$','FontUnits','points','interpreter','latex','FontSize',12,'FontName','Times')
ylabel('Amplication Ratio $X\frac{k}{F_{0}}$  ','FontUnits','points','interpreter','latex','FontSize',12,'FontName','Times')
% legend({'2'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Times','Location','NorthEast')
title('Amplitude','FontUnits','points','FontWeight','normal','FontSize',16,'FontName','Times')
grid on
print -depsc2 myplot.eps
print -dpdf myplot.pdf
hold off

