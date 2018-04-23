function FedV = FedV_2322(ee,lambda,G,L,H,xi,eta)
% The vector of elastic forces LdV based on St Venant Kircchoff material
% model
% Coded by MKM

II=diag([1,1]);
F=FF_2322(ee,L,H,xi,eta);
J=JJ_2322(ee,L,H,xi,eta);
dEEde=dEde_2322(ee,L,H,xi,eta);

% Cauchy-Green deformation tensor
CC=F'*F;
%Cinv=CC^(-1);

% Green-Lagrange strain
EE=zeros(2,2);
EE=1/2*(F'*F-II);

%2nd Piola Kirchhoff stress tensor
SS=zeros(2,2);
SS=lambda*II*trace(EE)+2*G*EE;   

% The vector of elastic forces without volume integration
FedV=zeros(1,12);
for kk=1:12,  
    FedV(kk)=0;
    for ii=1:2,
      for jj=1:2,           
            FedV(kk)=FedV(kk)+SS(ii,jj)*dEEde(ii,jj,kk);           
      end
   end
end
 

