function [ qb ] = AccelationSystem( y, model_def, alfa, beta )
%% this function is for accelation calculation
n = length(model_def.bodies);
M = Mass(model_def);
[ C, Cq, Cp, G ] = Dynamic_Constraint( y, model_def);
[ Qg ] = GravityForces( model_def );

[nn,~] = size(Cq);

LHS = [M, Cq'
      Cq, zeros(nn)];

RHS = [Qg
       G - 2*alfa*Cp - beta*beta*C];

x = LHS\RHS;

qb = x(1:3*n);

end