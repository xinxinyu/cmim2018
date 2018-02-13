function [ T, Y ] = odeRK4( fun, tspan, y0 )
%   Simple integration using Runge-Kutta Methods
%   fun - function handle with interface fun(t, y)
%   tspan - two element vector with dt and tend
%   y0 - initial conditions

% time step
dt = tspan(1);

T = 0:dt:tspan(2);

n = length(T);

% Y = zeros(length(y0), n);

% populate initial conditions
Y(:, 1) = y0(:);

% compute the solution
for i = 2:n    
    
    K1 = fun(T(i-1), Y(:, i-1));
    K2 = fun(T(i-1)+dt/2, Y(:, i-1)+0.5*dt*K1);
    K3 = fun(T(i-1)+dt/2, Y(:, i-1)+0.5*dt*K2);
    K4 = fun(T(i), Y(:, i-1)+dt*K3);
  
    Y(:, i) = Y(:, i-1) + dt * (K1+2*K2+2*K3+K4)/6;
end
end