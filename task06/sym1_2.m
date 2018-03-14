clc
clear all
syms alpha beta gama

%% Euler Rotation matrix
A1 = [cos(alpha) -sin(alpha) 0
      sin(alpha) cos(alpha) 0
      0 0 1];
A2  = [1 0 0
       0 cos(beta) -sin(beta)
       0 sin(beta) cos(beta)];
 
A3  = [cos(gama) -sin(gama) 0
       sin(gama) cos(gama) 0
       0 0 1];

A = A1*A2*A3;
AA1 = vpa(subs(A,[alpha beta gama],[45 45 45]));
AA2 = vpa(subs(A,[alpha beta gama],[90 30 -90]));

%% Euler parameters

e0=cos(beta/2)*cos((alpha+gama)/2);
e1=sin(beta/2)*cos((alpha-gama)/2);
e2=sin(beta/2)*sin((alpha-gama)/2);
e3=cos(beta/2)*sin((alpha+gama)/2);
e=[e0 e1 e2 e3];
ee = vpa(subs(e,[alpha beta gama],[90 30 -90]));




