function [ T, Y] = odeSemiFE( fun, tspan, y0)
%ODEFE Simple integration using Semi-implicit Euler Method
%   fun - function handle with interface fun(t, y)
%   tspan - two element vector with dt and tend
%   y0 - initial conditions

% time step
dt = tspan(1);

T = 0:dt:tspan(2);

n = length(T);

% Y = zeros(length(y0), n);

Y(:, 1) = y0(:);
V(:, 1) = y0(:);
for i = 2:n                                  
    
    Y(:, i) = Y(:, i-1) + dt * fun(T(i-1), V(:, i-1));
    V(:, i) = V(:, i-1) + dt * fun(T(i-1), Y(:, i));
end


end



