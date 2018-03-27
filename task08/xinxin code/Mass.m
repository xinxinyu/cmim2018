function [ MM ] = Mass( m, L, n );

MM = zeros(3*n);
for i = 1:n
    Ic(i) = 1/12 * m(i) * L(i) * L(i);
    M = diag([m(i), m(i), Ic(i)]);  
    MM(3*i-2:3*i,3*i-2:3*i) = M;
end

end
