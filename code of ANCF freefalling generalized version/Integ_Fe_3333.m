function Fe = Integ_Fe_3333(ElemDofs,E,nu,G,ks2,ks3,H,W,L,ee)
% Function: Integrate elastic forces using Gaussian quadrature

Fe0=zeros(1,ElemDofs);
Fev=zeros(1,ElemDofs);

% Get quadrature points and weights
nxi = 2;
neta = 2;
nzeta = 2;
 
[xiv,wxi]=gauleg2(-1,1,nxi);
[etav,weta]=gauleg2(-1,1,neta);
[zetav,wzeta]=gauleg2(-1,1,nzeta);
         
nxiL=2;

[xivL,wxiL] = gauleg2(-1,1,nxiL);

% Compute quadratures
for kk1=1:nzeta       
for jj1=1:neta
for ii1=1:nxi
    xxi=xiv(ii1)*L/2;
    eeta=etav(jj1)*H/2;
    zzeta=zetav(kk1)*W/2;
    Fe0=Fe0+Fe0dx_3333c(E,G,L,ks2,ks3,xxi,eeta,zzeta,ee)'*wxi(ii1)*weta(jj1)*wzeta(kk1);
end
end
end

for iiL=1:nxiL
    xxi=xivL(iiL)*L/2;
    Fev=Fev+FeVdx_3333c(E,nu,L,H,W,xxi,ee)'*wxiL(iiL);
end

% Sum up components (with scaling factors for interval change)
detF0=1/8*L*H*W;
detFv=L/2;
Fe=(Fe0*detF0+Fev*detFv);




