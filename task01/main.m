clc
close all
clear variables
tic
low_limit=0;
up_limit=pi;
no_splits=10e7;
fun=@(x)sin(x);

result = integral_trapezoid(...
    fun, low_limit, up_limit, no_splits )
toc
x=0:pi/100:pi;
y = fun(x);
trapz(x,y)