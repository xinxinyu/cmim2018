%% test mass spring system

m=1;
k=100;

dt = 1e-6;
tk = 10;

y0=[1 0];

[T, Y] = odeFE(@(t,y) Spring(t,y,m,k), [dt, tk], y0);