function [ C, Cq, Cp, G ] = Constr4bar( X,dX,t )


l1 = 0.2;
l2 = 0.4;
l3 = 0.3;
d1 = 0.35;
d2 = 0.1;
phinia1 = 2.36;
dphinia1 = 6.28;


theta1 = X(1);
theta2 = X(2);
theta3 = X(3);
dtheta1 = dX(1);
dtheta2 = dX(2);
dtheta3 = dX(3);

C = [l1*cos(theta1)+l2*cos(theta2)-l3*cos(theta3)-d1;
     l1*sin(theta1)+l2*sin(theta2)-l3*sin(theta3)-d2;
     theta1-phinia1-dphinia1*t];


Cq = [ -l1*sin(theta1), -l2*sin(theta2),  l3*sin(theta3);
        l1*cos(theta1),  l2*cos(theta2), -l3*cos(theta3);
                    1,              0,             0];

Ct = [0;
      0;
      -dphinia1];

Cp = Cq^(-1)*-Ct;

Ctt = [0;
       0;
       0];
   
Cqqq = [- l1*cos(theta1)*dtheta1^2 - l2*cos(theta2)*dtheta2^2 + l3*cos(theta3)*dtheta3^2;
        - l1*sin(theta1)*dtheta1^2 - l2*sin(theta2)*dtheta2^2 + l3*sin(theta3)*dtheta3^2; 
                                                                                       0];
Cqt = [0;
       0;
       0];
   
G = Cq^(-1)*(-Ctt-Cqqq-Cqt);
   
end