syms m1 m2 m3 m4 real
syms l1 l2 l3 l4 real
syms x1 y1 phi1 x2 y2 phi2 x3 y3 phi3 x4 y4 phi4 g real


body1.type = 'slender_link';
body1.m = m1;
body1.uix = -l1/2;
body1.uiy = 0;
body1.ujx = l1/2;
body1.ujy = 0;

body2.type = 'slender_link';
body2.m = m2;
body2.uix = -l2/2;
body2.uiy = 0;
body2.ujx = l2/2;
body2.ujy = 0;

body3.type = 'slender_link';
body3.m = m3;
body3.uix = -l3/2;
body3.uiy = 0;
body3.ujx = l3/2;
body3.ujy = 0;

body4.type = 'slender_link';
body4.m = 2;
body4.uix = -l4/2;
body4.uiy = 0;
body4.ujx = l4/2;
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

joint5.type = 'revolute';
joint5.ibody = 4;
joint5.jbody = 0;


model_def.joints = [joint1, joint2, joint3, joint4, joint5];

model_def.g = -g;

% initial condition of the system      
syms dx1 dy1 dphi1 dx2 dy2 dphi2 dx3 dy3 dphi3 dx4 dy4 dphi4 real
y0= [x1 y1 phi1 x2 y2 phi2 x3 y3 phi3 x4 y4 phi4];
y = [x1 y1 phi1 x2 y2 phi2 x3 y3 phi3 x4 y4 phi4 dx1 dy1 dphi1 dx2 dy2 dphi2 dx3 dy3 dphi3 dx4 dy4 dphi4]';
q = [x1 y1 phi1 x2 y2 phi2 x3 y3 phi3 x4 y4 phi4];
qdot= [dx1 dy1 dphi1 dx2 dy2 dphi2 dx3 dy3 dphi3 dx4 dy4 dphi4]';
[ C, Cq, Cp, G ] = Dynamic_Constraint( y, model_def);
CCq = jacobian (C,y0);

CM=Cq*qdot;
CC = jacobian(CM,q);
GG= -CC*qdot;



