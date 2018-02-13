
clear all
clc;
dt = 1e-5;
tk = 1;

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

eig(K,M)

y0=[0 0 3 0 0 0];


[T, Y1] = ode45(@(t,y) ThreeSpring(t,y,M,K), [0:dt:tk], y0);

[T, Y2] = ode15s(@(t,y) ThreeSpring(t,y,M,K), [0:dt:tk], y0);

[T, Y3] = odeRK4(@(t,y) ThreeSpring(t,y,M,K), [dt tk], y0);
Y3=Y3';

[T, Y4] = odeSemiFE(@(t,y) ThreeSpring(t,y,M,K), [dt tk], y0,K,M);
Y4=Y4';


plot(Y1(:,1))
hold on
plot(Y2(:,1))
hold on
plot(Y3(:,1))
hold on
plot(Y4(:,1))
hold off
