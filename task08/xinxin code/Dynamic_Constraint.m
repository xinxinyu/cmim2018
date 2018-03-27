function [ C, Cq, Cp, G ] = Dynamic_Constraint( y, L, n)

for i = 1:n
    xx(1,i) = y(3*i-2);
    yy(1,i) = y(3*i-1);
    phi(1,i) = y(3*i);
    dxx(1,i) = y(3*n+3*i-2);
    dyy(1,i) = y(3*n+3*i-1);
    dphi(1,i) = y(3*n+3*i);
end
%% Calculate the Constraint equation
C_revolute_direct1 = [xx(1)-L(1)/2*cos(phi(1))
                      yy(1)-L(1)/2*sin(phi(1))];
 for i = 1:n-1
     C3(i) = xx(1,i)+L(i)/2*cos(phi(i))-xx(1,i+1)+L(i+1)/2*cos(phi(i+1));
     C4(i) = yy(1,i)+L(i)/2*sin(phi(i))-yy(1,i+1)+L(i+1)/2*sin(phi(i+1));
 end    
C_revolute(1:2:2*n-3,:) = C3.';
C_revolute(2:2:2*n-2,:) = C4.';
% C_revolute_direct2 = [xx(n)+L(n)/2*cos(phi(n))-d(1)
%                       yy(n)+L(n)/2*sin(phi(n))-d(2)];
C = [C_revolute_direct1
     C_revolute];
%      C_revolute_direct2

%% Calculate the Jacobian equation
Cq_revolute_direct1 = [1 0 L(1)/2*sin(phi(1));
                       0 1 -L(1)/2*cos(phi(1))];
Cq(1:2,1:3) = Cq_revolute_direct1;                 

 for i=1:n-1
     Cq_revolute = [1 0 -L(i)/2*sin(phi(i)) -1 0 -L(i+1)/2*sin(phi(i+1))
                    0 1 L(i)/2*cos(phi(i)) 0 -1 L(i+1)/2*cos(phi(i+1))];   
  Cq(2*i+1:2*i+2,3*i-2:3*i+3) =  Cq_revolute;           
 end
% Cq_revolute_direct2 = [1 0 -L(n)/2*sin(phi(n));
%                        0 1 L(n)/2*cos(phi(n))];
% Cq(2*n+1:2*n+2,3*n-2:3*n) = Cq_revolute_direct2;

%% Calculate Cp
qp = y(3*n+1:6*n);
Cp = Cq*qp;

%% Calculate G

G_revolute_direct1 = [-(L(1)*dphi(1)^2*cos(phi(1)))/2
                      -(L(1)*dphi(1)^2*sin(phi(1)))/2];
 for i = 1:n-1
     G3(i) = (L(i)*dphi(i)^2*cos(phi(i)))/2+(L(i+1)*dphi(i+1)^2*cos(phi(i+1)))/2;
     G4(i) = (L(i)*dphi(i)^2*sin(phi(i)))/2+(L(i+1)*dphi(i+1)^2*sin(phi(i+1)))/2;
 end    
G_revolute(1:2:2*n-3,:) = G3.';
G_revolute(2:2:2*n-2,:) = G4.';
% G_revolute_direct2 = [(L(n)*dphi(n)^2*cos(phi(n)))/2
%                       (L(n)*dphi(n)^2*sin(phi(n)))/2];
G = [G_revolute_direct1
     G_revolute];
%      G_revolute_direct2];

end
