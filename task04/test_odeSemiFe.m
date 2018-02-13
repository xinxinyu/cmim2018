%% test mass spring system
clc 
clear
close all

m1=2;
m2=1.5;
m3=1;
m=diag([m1 m2 m3]);
k1=20000;
k2=15000;
k3=10000;
k=[k1+k2 -k2 0;
    -k2 k2+k3 -k3;
      0 -k3 k3];
    

dt = 1e-5;
tk = 10;

y0=[0 0 0.1 0 0 0];

[T, Y] = odeSemiFE(@(t,y) Spring(t,y,m,k), [dt, tk], y0,k,m);

