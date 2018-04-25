function xlocAll= xlocAllANCF_2322(nloc)
% function xlocAllANCF makes full xloc for all elements 
% row - element, column - dof


[nl,m] = size(nloc);

xlocAll = zeros(nl,4*3);

for k = 1:nl
  n1 = [nloc(k,1)*4-3:1:nloc(k,1)*4];
  n2 = [nloc(k,2)*4-3:1:nloc(k,2)*4];
  n3 = [nloc(k,3)*4-3:1:nloc(k,3)*4];
  xlocAll(k,:) = [n1 n2 n3];
end