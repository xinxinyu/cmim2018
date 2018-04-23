function [ C ] = ConstraintEquation( x1, y1, phi1, x2, y2, phi2, u1_1, u1_2, u2_1, u2_2)

   
C1 = x1+u1_1*cos(phi1)-u1_2*sin(phi1)-x2-u2_1*cos(phi2)+u2_2*sin(phi2);
C2 = y1+u1_1*sin(phi1)+u1_2*cos(phi1)-y2-u2_1*sin(phi2)-u2_2*cos(phi2);
C = [C1 C2]';

end




