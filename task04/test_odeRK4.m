%% test if function works

% test function  y'=r*y

m=1;
k=100;

dt = 1e-4;
tk = 10;

y0=[1 0];


[T, Y] = ode45(@(t,y) Spring(t,y,m,k), [0:dt:tk], y0);
% [T, Y] = odeRK4(@(t,y) Spring(t,y,m,k), [dt, tk], y0);

plot(Y);