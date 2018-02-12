%% test mass spring system
clc 
clear
close all

m=1;
k=100;

dt = 1e-3;
tk = 10;

y0=[1 0];

[T, Y] = odeSemiFE(@(t,y) Spring(t,y,m,k), [dt, tk], y0);

