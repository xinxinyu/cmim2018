function C = ConsEqua4Bar(X,t)

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

C = [l1*cos(theta1)+l2*cos(theta2)-l3*cos(theta3)-d1;
     l1*sin(theta1)+l2*sin(theta2)-l3*sin(theta3)-d2;
     theta1-phinia1-dphinia1*t];

end