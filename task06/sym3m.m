clc 
clear all
syms x y z phi1 phi2 phi3 t
syms dx dy dz dphi1 dphi2 dphi3 
syms ddx ddy ddz ddphi1 ddphi2 ddphi3 

q = [x;y;z;phi1;phi2;phi3];
dq = [dx;dy;dz;dphi1;dphi2;dphi3];
ddq = [ddx;ddy;ddz;ddphi1;ddphi2;ddphi3];
C1 = [x^2+y+sqrt(z)+sin(phi1)];
C2 = [x*y+x*z+y*sin(phi3)+t^3];
C3 = [sin(phi2)+x^1.5+t];
C = [C1 
     C2 
     C3];
Cq = jacobian (C,q);
Ct = diff(C,t);
dC = Cq*dq+Ct;

Cqdq = jacobian(Cq*dq,q);
Cqt = diff(Cq ,t);
Ctt = diff(C,t,2);
ddC = Cq*ddq+Cqdq*dq+2*Cqt*dq+Ctt;                