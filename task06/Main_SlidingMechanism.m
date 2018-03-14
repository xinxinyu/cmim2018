clc
clear all

% syms r phi h0 c real

c = 3;
h0 = 4;

fun = @(phi) c/cos(phi)-h0/sin(phi);
dfun = @(phi) c*sin(phi)/(cos(phi))^2+h0*cos(phi)/(sin(phi))^2;
phi0 = pi/4;
[xr, flag] = NewtonRaphson(fun, dfun, phi0);
r= c/cos(xr)