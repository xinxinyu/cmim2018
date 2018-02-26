
syms x1 x2 y1 y2


fun=@(x)x^2;

x = [x1 x2];
y1 = fun(x1);
y2 = fun(x2);

interval = [0 1000];

func = [y1 y2];

X = BruteForce(func, x1 , x2 ,interval);