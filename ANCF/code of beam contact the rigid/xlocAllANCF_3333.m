function xlocAll= xlocAllANCF_3333(nloc)
% function xlocAllANCF makes full xloc for all elements 
% row - element, column - dof


[nl,m] = size(nloc);

xlocAll = zeros(nl,9*3);

for k = 1:nl
  n1 = [nloc(k,1)*9-8:1:nloc(k,1)*9];
  n2 = [nloc(k,2)*9-8:1:nloc(k,2)*9];
  n3 = [nloc(k,3)*9-8:1:nloc(k,3)*9];
  xlocAll(k,:) = [n1 n2 n3];
end