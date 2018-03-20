function Cq = Jaco4Bar(X,t)



l1 = 0.2;
l2 = 0.4;
l3 = 0.3;
theta1 = X(1);
theta2 = X(2);
theta3 = X(3);

Cq = [ -l1*sin(theta1), -l2*sin(theta2),  l3*sin(theta3);
        l1*cos(theta1),  l2*cos(theta2), -l3*cos(theta3);
                    1,              0,             0];
end