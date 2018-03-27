clear all
close all
clc

%% Define mass matrix

%Define Mass
% n = 3;
% m1 = 2;
% m2 = 4;
% m3 = 3;
% m = [m1 m2 m3];
% l1 = 0.2;
% l2 = 0.4;
% l3 = 0.3;
% d1 = 0.35;
% d2 = 0.1;
% d = [d1 d2];
% body1.coord = [-0.07,0.07,2.36];
% body2.coord = [0.03,0.25,0.57];
% body3.coord = [0.275,0.07,0.23];
n = 4;
m1 = 2;
m2 = 4;
m3 = 3;
m4 = 2;
m = [m1 m2 m3 m4];
l1 = 0.2;
l2 = 0.4;
l3 = 0.3;
l4 = 0.2;
l = [l1 l2 l3 l4];
g = -9.81;
M = Mass(m,l,n);
Qg = GravityForces( m, g, n );

%% Solve the system
odefun = @(t, y) [y(3*n+1:6*n); M\Qg];

body1.coord = [l1/2,0,0];
body2.coord = [l1+l2/2,0,0];
body3.coord = [l1+l2+l3/2,0,0];
body4.coord = [l1+l2+l3+l4/2,0,0];

yy0 = [body1.coord body2.coord body3.coord body4.coord];
dy0 = zeros (1,length(yy0));
y0 = [yy0 dy0]';

tspan = linspace(0, 1, 1001);
[T, Y] = ode45(odefun, tspan, y0);

%% Symbolic computation
sy = sym('y', [6*n, 1]);
sL = sym('l', [n, 1]);
% sd = sym('d', [2, 1]);
[ C, Cq, Cp, G ] = Dynamic_Constraint( sy, sL, n);

%% Verify G
q = sy(1:3*n);
qdot = sy(3*n+1:6*n);
CM=Cq*qdot;
CC = jacobian(CM,q);
GG= -CC*qdot;
error=GG-G;

%% Solve MBD system


alfa = 1;
beta = 1;

odefun2 = @(t, y) [y(3*n+1:6*n)
    AccelationSystem( y, M, l, g, alfa, beta, n )];

[T, YY] = ode45(odefun2, tspan, y0);
% 
% plot(T, Y(:, 2), T, Y(:, 5));


% m1 = 2;
% l1 = 0.2;
% M1 = MassMatrix(m1, l1);
% 
% m2 = 4;
% l2 = 0.4;
% M2 = MassMatrix(m2, l2);
% 
% m3 = 3;
% l3 = 0.3;
% M3 = MassMatrix(m3, l3);

% 
% 
% 
% 
% 
% n = 3;
% body1.coord = [-0.07,0.07,2.36];
% body2.coord = [0.03,0.25,0.57];
% body3.coord = [0.275,0.07,0.23];
% % body4.coord = [3,4,5];

% % l4 = 1;
% d1 = 0.35;
% d2 = 0.1;
% 
% L = [l1 l2 l3];
% y0 = [body1.coord body2.coord body3.coord];
% dy0 = zeros (1,length(y0));
% y = [y0 dy0];
% T=0; 










% 
% 
% % Defining the generalized coordinates
% syms R11 R21 theta1 R12 R22 theta2 R13 R23 theta3
% syms dR11 dR21 dtheta1 dR12 dR22 dtheta2 dR13 dR23 dtheta3 
% syms m1 m2 l1 l2 m3 b h  Torque t
% syms ctheta3 k g R13ini
% 
% q1 = [R11;R21;theta1];
% q2 = [R12;R22;theta2];
% q3 = [R13;R23;theta3];
% q = [q1;q2;q3];
% 
% dq1 = [dR11;dR21;dtheta1];
% dq2 = [dR12;dR22;dtheta2];
% dq3 = [dR13;dR23;dtheta3];
% dq = [dq1;dq2;dq3];
% 
% % Rotation matrix of the bodies
% 
% A1 = [cos(theta1),-sin(theta1);
%       sin(theta1), cos(theta1)];
%   
% A2 = [cos(theta2),-sin(theta2);
%       sin(theta2), cos(theta2)];
% 
% A3 = [cos(theta3),-sin(theta3);
%       sin(theta3), cos(theta3)]; 
%   
%    
% % Forming the system constraint equations.
% 
% X0 = [0;0];  % Ground - crank revolute joint attacehment point.
% 
% C(1:2) = [R11;R21] + A1*[-l1/2;0] - X0; % Ground - crank revolute joint
% 
% C(3:4) = [R11;R21] + A1*[l1/2;0] - ([R12;R22] + A2*[l2/2;0]); % Crank-shaft revolute joint.
% 
% C(5:6) = [R12;R22] + A2*[-l2/2;0] - ([R13;R23] + A3*[0;0]); % Shaft - block revolute joint.
% 
% C(7) = theta3 - ctheta3;
% 
% hunit = [0;1]; % Vector perpendicular to the joint axis.
% rp = [R13;R23];
% 
% C(8) = hunit.'*rp;
% 
% 
% 
% C = C(:);
% 
% % Forming the jacobian matrix
% 
% Cq = jacobian(C,q);
% 
% %% Forming the mass matrices of the bodies
% 
% I1 = 1/12*m1*l1^2;
% I2 = 1/12*m2*l2^2;
% I3 = 1/12*m3*(b^2+h^2);
% 
% 
% M1 = [m1,  0,  0;
%        0, m1,  0;
%        0,  0, I1];
%    
% M2 = [m2,  0,  0;
%        0, m2,  0;
%        0,  0, I2];
% 
% M3 =  [m3,  0,  0;
%        0, m3,  0;
%        0,  0, I3];
%    
% M(1:3,1:3) = M1;
% M(4:6,4:6) = M2;
% M(7:9,7:9) = M3;
% 
% %% Forming the system matrix
% 
% SysM = [M ,                            Cq.';
%         Cq, zeros(size(Cq,1),size(Cq.',2))];
%     
% % converting the system matrix into a function file.
%     
%     
% %matlabFunction(SysM,'file','SysMn');
% 
% 
% %% Vector of generalized external forces
% 
% Qe1 = [0;-m1*g;Torque];
% Qe2 = [0;-m2*g;0];
% Qe3 = [-k*(R13-R13ini);-m3*g;0];
% 
% Qe = [Qe1;Qe2;Qe3];
% 
% % Vector of constraint forces
% 
% Qd = -diff(C,t,2) - jacobian(Cq*dq,q)*dq - 2*diff(Cq,t)*dq;
% 
% 
% %% Forming the system vector
% 
% SysF = [Qe;Qd];
% 
% %matlabFunction(SysF,'file','SysFn');
% 
% %% EOM without stabilization
% 
% 
% EOM = simplify(SysM\SysF);
% 
% %% EOM with stabilization
% 
% 
% %Baumgarten stabilisation applied to Langrange multipliers NB! Maybe not the best parameters
% deltat=0.01;   % 1st parameter of Baumgarte
% alpha=5000;    % 2nd parameter of Baumgarte
% beta=5000;     % 3rd parameter of Baumgarte
% Qd=Qd-(2*alpha*Cq*dq + beta^2*C);
% 
% SysF = [Qe;Qd];    
% 
% EOMs = simplify(SysM\SysF);
% 
% matlabFunction(EOM,'file','EOMn');
% matlabFunction(EOMs,'file','EOMsn');
% 
% %% Integration of the equation of motion
% 
% % Vector of initial conditions
% 
% theta1ini = 30*pi/180;
% 
% l1ini = 0.1;
% l2ini = 0.3;
% 
% theta2ini = pi - asin((sin(theta1ini)*l1ini/l2ini));
% theta2inid = theta2ini*180/pi; % Just to validate the calculation.
% 
% R1ini = double(subs(A1*[l1/2;0],[theta1,l1],[theta1ini,l1ini]));
% R2ini = 2*R1ini + double(subs((A2*[-l2/2;0]),[theta2,l2],[theta2ini,l2ini]));
% R3ini = 2*R1ini + double(subs((A2*[-l2;0]),[theta2,l2],[theta2ini,l2ini]));
% 
% y0 = [R1ini;theta1ini;R2ini;theta2ini;R3ini;0;zeros(9,1)];
% % List of system parameters
% 
% par(1) = 3;                 % m1
% par(2) = 3;                 % m2
% par(3) = 0.5;               % m3
% par(4) = l1ini;             % l1
% par(5) = l2ini;             % l2
% par(6) = 9.81;              % g
% par(7) = 0.01;              % h               
% par(8) = 0.01;              % b
% par(9) = 0;                 % ctheta3;
% par(10) = 1000;              % k   
% par(11) = R3ini(1);         % R13ini
% 
% par = par(:);
% 
% 
% tini = 0;
% dt = 0.01;
% tend = 5;
% 
% tspan = tini:dt:tend;
% 
% options = [];
% 
% % Use stabilization in the EOM
% 
% %stabilization = 'yes';
% stabilization = 'no';
% 
% [t,y] = ode45(@eomDDrivenCS, tspan, y0,options,par,stabilization );
