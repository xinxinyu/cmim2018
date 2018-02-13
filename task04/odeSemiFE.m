function [ T, Y] = odeSemiFE( fun, tspan, y0,K,M)
%   Simple integration using Semi-implicit Euler Method
%   fun - function handle with interface fun(t, y)
%   tspan - two element vector with dt and tend
%   y0 - initial conditions

% time step
dt = tspan(1);

T = 0:dt:tspan(2);

n = length(T);

% Y = zeros(length(y0), n);

Y(:, 1) = y0(:);

V(1:3, 1) = y0(4:6)';
V(4:6, 1) = -1*(M^-1)*K*y0(1:3)';

for i = 2:n                                  
    
    V(:, i) = V(:, i-1) + dt * fun(T(i-1), V(:, i-1));
    Y(:, i) = Y(:, i-1) + dt * V(:, i);
end


end



