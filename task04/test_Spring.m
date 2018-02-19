%% test mass spring system

m=1;
k=100;

dt = 1e-7;
tk = 1;

y0=[1 0];

[T, Y] = odeFE(@(t,y) Spring(t,y,m,k), [dt, tk], y0);


V=0.5*k*Y(1,:).^2;
U=0.5*m*Y(2,:).^2;
T=V+U;