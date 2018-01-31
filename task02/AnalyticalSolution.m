syms m k c

syms x(t)

xp = diff(x, t);
xpp = diff(xp, t);

Eq = m*xpp + k*x+c*xp;

solution_free = dsolve(Eq == 0);

disp(solution_free)