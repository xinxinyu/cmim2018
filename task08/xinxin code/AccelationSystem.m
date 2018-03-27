function [ qb ] = AccelationSystem( y, M, l, g, alfa, beta, n )

M = Mass(M, l, n);
[ C, Cq, Cp, G ] = Dynamic_Constraint( y, l, n);
[ Qg ] = GravityForces( M, g, n );

[nn,~] = size(Cq);

LHS = [M, Cq'
      Cq, zeros(nn)];

RHS = [Qg
       G - 2*alfa*Cp - beta*beta*C];

x = LHS\RHS;

qb = x(1:3*n);

end