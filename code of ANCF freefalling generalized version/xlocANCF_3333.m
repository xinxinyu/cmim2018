function loc = xlocANCF_3333(nodes,comps)
% make location vector in x from nodelist and componentlist

loc = [];
for n=1:length(nodes)
  nn = nodes(n);
  loc = [loc 9*nn-9+comps];
end
  