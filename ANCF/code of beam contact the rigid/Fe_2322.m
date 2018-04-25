function Fe = Fe_2322(ElemDofs,lambda,G,L,H,W,ee)
% Function: Integrate elastic forces using Gaussian quadratures
% Coded by MKM
Fe=zeros(1,ElemDofs);

detF0=L*H/4;

% number of integration points
nxi=3;
neta=3;

% [xiv,wxi]=gauleg2(-1,1,nxi);
% [etav,weta]=gauleg2(-1,1,neta);

% Ip:s and weights are given here instead of use of gauleg2 functiono 
% to compute them for every single element. Should increase computational 
% efficiency

xiv=[-0.774596669241483, 0, 0.774596669241483];
etav=xiv;

wxi=[0.555555555555555, 0.888888888888889, 0.555555555555555];
weta=wxi;

for jj1=1:neta,        
    for ii1=1:nxi,
        xi=xiv(ii1);
        eta=etav(jj1);
        Fe=Fe+FedV_2322(ee,lambda,G,L,H,xi,eta)*wxi(ii1)*weta(jj1);   
    end
end

Fe=Fe*detF0*W;



