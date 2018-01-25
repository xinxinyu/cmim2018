clc
close all
clear variables

low_limit=0;
up_limit=pi;
no_splits=100;
fun=@(x)sin(x);

result = integral_trapezoid(...
    fun, low_limit, up_limit, no_splits )

x=0:pi/100:pi;
y = fun(x);
trapz(x,y)