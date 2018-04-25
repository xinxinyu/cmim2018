function [ T, u, v] = odeSemiFE_constant_contactforce(tspan, u0, v0,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse,mu)
% nx   DOFs 
% ndof DOFs after linear constraints
T = tspan;
n = length(T);
u = zeros(length(u0),n);
v = zeros(length(v0),n);
u(:,1) = u0;
v(:,1) = v0;

for l = 1:n-1
    dt = T(l+1)-T(l);                                      %% time step
     
    ee=zeros(nx,1);         % positions, systeemin kaikki koordinaatit, sek� lineaarisesti rajoitetut
    eedot=zeros(nx,1);      % velocities, -"-

% positions
ee(bc)=u(1:ndof);
% velocities
eedot(bc)=v(1:ndof); 

if bcInd~=0
    ee(bcInd)=ee0bc;
    eedot(bcInd)=ee0dotbc;
end

ee=ee(:);
eedot=eedot(:);

%% Calculation about external and elastic forces
Fe=zeros(nx,1);
Fext=zeros(nx,1);

% loop over the elements - calculation of elastic and external forces
% (e.g. gravity)
for k = 1:nl
        xlock = xloc(k,:);
        eek=ee(xlock);
      
        Fek=Fe_2322(ElemDofs,lambda,G,L,H,W,eek)'; 
        Fextk=FextBF(rho,g,L,H,W);
        
        Fe(xlock)=Fe(xlock)+Fek;
        Fext(xlock)=Fext(xlock)+Fextk;    
end

% Elimination of linear constraints
Fextc=Fext(bc);
Fec=Fe(bc);

%% Calculation about spring force and damping force


    SP1 = [1 0 -H/2 0 0 0 0 0 0 0 0 0 ;
           0 1 0 -H/2 0 0 0 0 0 0 0 0 ]; 
    SP2 = [0 0 0 0 1 0 -H/2 0 0 0 0 0 ;
           0 0 0 0 0 1 0 -H/2 0 0 0 0 ];                    %% shape function at the contact point P     %% Xinxin rectified it
    SP3 = [0 0 0 0 0 0 0 0 1 0 -H/2 0 ;
           0 0 0 0 0 0 0 0 0 1 0 -H/2 ];

 rP1(:,l) = SP1*u(:,l);                                   %% position value at the contact point P     %% Xinxin rectified it
 Phi1 = rP1(2,l);                                         %% gap function for CCP problem              %% Xinxin rectified it       
 rP2(:,l) = SP2*u(:,l); 
 Phi2 = rP2(2,l);
 rP3(:,l) = SP3*u(:,l);                                   %% position value at the contact point P     %% Xinxin rectified it
 Phi3 = rP3(2,l); 

c = 1500;
ContactSpring1 = CalculateForce(Phi1,c);
Spring1 = [0 ContactSpring1 0 0 0 0 0 0 0 0 0 0]';
ContactSpring2 = CalculateForce(Phi2,c);
Spring2 = [0 0 0 0 0 ContactSpring2 0 0 0 0 0 0]';
ContactSpring3 = CalculateForce(Phi3,c);
Spring3 = [0 0 0 0 0 0 0 0 0 ContactSpring3 0 0]';

v(:,l+1) = v(:,l) + (MInvSparse*Fextc-MInvSparse*Fec);   % Eq 13
%v(:,l+1) = v(:,l) + (MInvSparse*Fextc-MInvSparse*Fec + MInvSparse*(Spring1+Spring2+Spring3));   % Eq 13
u(:,l+1) = u(:,l) + v(:,l+1)*dt;          % Eq 9a


Fextdata(l,:)=Fext;
   
end
save Fextdata.mat Fextdata
