%% test if function works

% test function  y'=r*y

m=1;
k=100;

dt = 1e-3;
tk = 10;

y0=[1 0];

[T, Y] = odeRK4(@(t,y) Spring(t,y,m,k), [dt, tk], y0);
