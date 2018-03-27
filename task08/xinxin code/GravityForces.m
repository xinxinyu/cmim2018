function [ Qg ] = GravityForces( m, g, n )


for i = 1:n
    Q = [0
        m(i)*g
             0];  
    Qg(3*i-2:3*i,1) = Q;
end


end