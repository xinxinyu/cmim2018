function [ Qg ] = GravityForces( model_def )
%% this funciton is made to calculate the gravity force matrix
n = length(model_def.bodies);
for i = 1:n
    Q = [0
        model_def.bodies(i).m*model_def.g
         0];  
    Qg(3*i-2:3*i,1) = Q;
end


end