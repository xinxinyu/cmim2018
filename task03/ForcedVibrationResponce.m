clear all


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
plot(A1,'LineWidth', 1.5)
hold on
plot(A2,'LineWidth', 1.5)
hold on
plot(A3,'LineWidth', 1.5)
hold on
plot(A4,'LineWidth', 1.5)
hold on
plot(A5,'LineWidth', 1.5)
hold on
plot(A6,'LineWidth', 1.5)
hold on
plot(A7,'LineWidth', 1.5)
xlim([0,100]);
ylim([0,5]);
hold off
xlabel('Frequency Ratio $$r={f\over f}$$','FontUnits','points','interpreter','latex','FontSize',12,'FontName','Times')
ylabel('Amplication Ratio $$X,{k\over F}$$','FontUnits','points','FontWeight','normal','FontSize',12,'FontName','Times')
% legend({'2'},'FontUnits','points','interpreter','latex','FontSize',12,'FontName','Times','Location','NorthEast')
title('Amplitude','FontUnits','points','FontWeight','normal','FontSize',12,'FontName','Times')
grid on
print -depsc2 myplot.eps
print -dpdf myplot.pdf

hold off



