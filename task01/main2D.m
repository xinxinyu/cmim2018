clc
close all
clear variables

tic
low_limitX=-3;
up_limitX=3;
low_limitY=-5;
up_limitY=5;
no_splits=100;


 result = integral_trapezoid(@(y)(integral_trapezoid(@(x)(-x^4-y^4),low_limitX, up_limitX,no_splits)),low_limitY, up_limitY,no_splits);

toc
x=-3:6/100:3;
y=-5:10/100:5;
[X,Y]=meshgrid(x,y);
fun=-X.^4-Y.^4;
I = trapz(y,trapz(x,fun,2));
surf(X,Y,fun)