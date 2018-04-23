function [ MM ] = Mass( model_def );
%% this function is made to calculate the general mass matrix
n = length(model_def.bodies);
MM = zeros(3*n);
for i = 1:n
    Ic(i) = 1/12 * model_def.bodies(i).m * model_def.bodies(i).l * model_def.bodies(i).l;
    M = diag([model_def.bodies(i).m, model_def.bodies(i).m, Ic(i)]);  
    MM(3*i-2:3*i,3*i-2:3*i) = M;
end

end
