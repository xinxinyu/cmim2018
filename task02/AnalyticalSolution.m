%% Symbolic toolbox tutorial - in appliction to SDOF systems
% Open help for symbolic toolbox.

%% use analytical solution

doc symbolic
%% Define symbolic variables
%%
syms m k c real
assumeAlso(m > 0);
assumeAlso(k > 0);
assumeAlso(c > 0);

syms x(t) position gap1
%% Analyze free vibrations
% $$m\ddot{x}+kx=0$$
%%
xp = diff(x, t);
xpp = diff(xp, t);

Eq = m*xpp + k*x;

solution_free = dsolve(Eq == 0);
c 
disp(solution_free)
%% rewrite
%%
oms = sqrt(k/m);
syms omega real

solution_free = subs(solution_free, oms, omega);

disp(solution_free)
%% Solve with initial conditions
% $$\left\{ \begin{array}{c}m\ddot{x}+kx=0\\x\left(0\right)=A,\,\dot{x}\left(0\right)=0\end{array}\right.$$
%%
syms A real
solution_free_ic = dsolve(Eq == 0, [xp(0) == 0, x(0) == A] );
solution_free_ic = subs(solution_free_ic, oms, omega);
disp(solution_free_ic)
% display equation variables
symvar(solution_free_ic)
% %% Plot solution
% %%
% figure
% fplot(subs(solution_free_ic, [A, omega], [-1, 2]), [0, 10])
% grid on

%% Get numerical solution
tsol = linspace(0, 10, 501);
xsol = subs(solution_free_ic, [A, omega, t], {-1, sqrt(10/1), tsol});
vsol = subs(diff(solution_free_ic, t), [A, omega, t], {-1, sqrt(10/1), tsol});

position=vpa(xsol)'; %% Get numerical result of the position


%% use ode function

format long
m=1;
k=10;
c=0;
A=100;
omga=2;

v0=[-1 0]; %% Initial position and velocity of the system
dt = 0.02;
TT=10;

%% use ode solver with reducing the errors
RelTol=1e-6;
AbsTol=1e-9; 
options = odeset('AbsTol',AbsTol,'RelTol',RelTol);
[t1,v1]=ode45(@(t,v) SysPosition(t,v,m,k,c,A,omga),[0 : dt : TT], v0,options);
[t2,v2]=ode15s(@(t,v) SysPosition(t,v,m,k,c,A,omga),[0 : dt : TT], v0,options);
[t3,v3]=ode23(@(t,v) SysPosition(t,v,m,k,c,A,omga),[0 : dt : TT], v0,options);

%% use ode solver without reducing the errors

[t4,v4]=ode45(@(t,v) SysPosition(t,v,m,k,c,A,omga),[0 : dt : TT], v0);
[t5,v5]=ode15s(@(t,v) SysPosition(t,v,m,k,c,A,omga),[0 : dt : TT], v0);
[t6,v6]=ode23(@(t,v) SysPosition(t,v,m,k,c,A,omga),[0 : dt : TT], v0);

%% Difference of the position between the analytical solution and ode solvers
gap1=position-v1(:,1);  % ode45 with reducing error
gap2=position-v2(:,1);  % ode15s with reducing error
gap3=position-v3(:,1);  % ode23 with reducing error
gap4=position-v4(:,1);  % ode45 without reducing error


%% Plot all the figures
% 
figure(1)   %% Plot the difference of the position between the analytical solution and ode solvers
plot(gap1,'LineWidth', 1.5 );
hold on
plot(gap2,'--','LineWidth', 1.5 );
hold on
legend('with reducing the errors','without reducing the errors');
plot(gap3,':','LineWidth', 1.5 );
grid on
legend('ode45','ode15s','ode23')
hold off

figure(2)   %% Plot the position and velocity with different ode solvers
plot(v1(:,2),'LineWidth', 1.5 );
hold on
plot(v1(:,1),'LineWidth', 1.5 );
hold on
plot(v2(:,2),'--','LineWidth', 1.5 );
hold on
plot(v2(:,1),'--','LineWidth', 1.5 );
hold on
plot(v3(:,2),':','LineWidth', 1.5 );
hold on
plot(v3(:,1),'--','LineWidth', 1.5 );
hold off
grid on
legend('velocity-ode45','position-ode45','velocity-ode15s','position-ode15s','velocity-ode23','position-ode23')

figure(3)   %%
plot(gap1,'LineWidth', 1.5 );
hold on
plot(gap4,'--','LineWidth', 1.5 );
hold on
legend('with reducing the errors','without reducing the errors');
grid on

























