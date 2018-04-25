function loc = xlocANCF_2322(nodes,comps)
% make location vector in x from nodelist and componentlist

loc = [];
for n=1:length(nodes)
  nn = nodes(n);
  loc = [loc 4*nn-4+comps];
end
  