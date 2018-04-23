%clear all
%close all
%clc
%% initial data
% model definition
syms m1 m2 m3 m4 real
syms l1 l2 l3 l4 real
syms x1 y1 phi1 x2 y2 phi2 x3 y3 phi3 x4 y4 phi4 real

body1.type = 'slender_link';
body1.m = 2;
body1.l = 2;
body1.uix = -1;
body1.uiy = 0;
body1.ujx = 1;
body1.ujy = 0;

body2.type = 'slender_link';
body2.m = 4;
body2.l = 4;
body2.uix = -2;
body2.uiy = 0;
body2.ujx = 2;
body2.ujy = 0;

body3.type = 'slender_link';
body3.m = 3;
body3.l = 3;
body3.uix = -1.5;
body3.uiy = 0;
body3.ujx = 1.5;
body3.ujy = 0;

body4.type = 'slender_link';
body4.m = 2;
body4.l = 2;
body4.uix = -1;
body4.uiy = 0;
body4.ujx = 1;
body4.ujy = 0;

model_def.bodies = [body1, body2, body3, body4];

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
joint4.jbody = 4;

model_def.joints = [joint1, joint2, joint3, joint4];

model_def.g = -9.81;

%% initial condition of the system
yy0 = [1 0 0 4 0 0 7.5 0 0 10 0 0];                        % initial position

dy0 = zeros (1,length(yy0));                               % initial velocity
y0 = [yy0 dy0]';

[ C, Cq, Cp, G ] = Dynamic_Constraint( y0, model_def);

% % simulation settings
% settings.method.type = 'baumgarte';
% settings.method.alfa = 1;
% settings.method.beta = 1;

%% Solve MBD system

y0 = [yy0 dy0]';
alfa = 1;
beta = 1;
tspan = linspace(0, 1, 1001);
n = length(model_def.bodies);
odefun2 = @(t, y) [y(3*n+1:6*n)
    AccelationSystem( y, model_def, alfa, beta)];

[T, YY] = ode45(odefun2, tspan, y0);










