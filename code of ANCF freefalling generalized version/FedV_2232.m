function FedV = FedV_2232(ee,L,H,W,E,G,lambda,K,xi,eta,zeta)
% The vector of elastic forces LdV based on St Venant Kircchoff material
% model
% Coded by MKM

II=diag([1,1]);
F=FF_2232(ee,L,H,W,xi,eta,zeta);
J=JJ_2232(ee,L,H,W,xi,eta,zeta);

dEEde=dEde_2232(ee,L,H,W,xi,eta,zeta);

% Cauchy-Green deformation tensor
CC=F'*F;

%Cinv=CC^(-1);

% Green-Lagrange strain
EE=zeros(2,2);
EE=1/2*(F'*F-II);

%2nd Piola Kirchhoff stress tensor
SS=zeros(2,2);
SS=lambda*II*trace(EE)+2*G*EE;   


FedV=zeros(1,12);

for kk=1:12,  
    FedV(kk)=0;
    for ii=1:2,
      for jj=1:2,           
            FedV(kk)=FedV(kk)+SS(ii,jj)*dEEde(ii,jj,kk);           
      end
   end
end
 

