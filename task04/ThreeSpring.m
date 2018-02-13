function dy = ThreeSpring(t,y)

m1=2;
m2=1.5;
m3=1;
M=diag([m1,m2,m3]);
k1=20000;
k2=15000;
k3=10000;
K=[k1+k2 -k2 0;
    -k2 k2+k3 -k3;
    0 -k3 k3];
Dist=[y(1) y(2) y(3)]';
M_1 = inv(M);

Accela=-M_1*K*Dist;

dy(1)=y(4);
dy(2)=y(5);
dy(3)=y(6);
dy(4)=Accela(1);
dy(5)=Accela(2);
dy(6)=Accela(3);

end