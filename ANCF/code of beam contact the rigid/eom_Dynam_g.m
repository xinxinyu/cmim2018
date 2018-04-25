function yddot = eom_Dynam_g(t,u,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse,D,di0,dt)
% nx   DOFs 
% ndof DOFs after linear constraints
%% CALCULATE the velocity of the element

ee=zeros(nx,1);         % positions, systeemin kaikki koordinaatit, sek� lineaarisesti rajoitetut

% positions
ee(bc)=u(1:ndof);

if bcInd~=0
    ee(bcInd)=ee0bc;
end

ee=ee(:);

Fe=zeros(nx,1);
Fext=zeros(nx,1);

%% loop over the elements - calculation of elastic and external forces
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
Fextc=Fext(bc);                     % external forces
Fec=Fe(bc);                         % elastic forces

    SP1 = [1 0 -H/2 0 0 0 0 0 0 0 0 0 ;
           0 1 0 -H/2 0 0 0 0 0 0 0 0 ]; 
    SP2 = [0 0 0 0 1 0 -H/2 0 0 0 0 0 ;
           0 0 0 0 0 1 0 -H/2 0 0 0 0 ];                    %% shape function at the contact point P     %% Xinxin rectified it
    SP3 = [0 0 0 0 0 0 0 0 1 0 -H/2 0 ;
           0 0 0 0 0 0 0 0 0 1 0 -H/2 ];

 rP1 = SP1*u;                                              %% position value at the contact point P     %% Xinxin rectified it
 Phi1 = rP1(2);                                             %% gap function for CCP problem              %% Xinxin rectified it       
 rP2 = SP2*u; 
 Phi2 = rP2(2);
 rP3 = SP3*u;                                              %% position value at the contact point P     %% Xinxin rectified it
 Phi3 = rP3(2); 

c = 150000;
ContactSpring1 = CalculateForce(Phi1,c);
Spring1 = [0 ContactSpring1 0 0 0 0 0 0 0 0 0 0]';
ContactSpring2 = CalculateForce(Phi2,c);
Spring2 = [0 0 0 0 0 ContactSpring2 0 0 0 0 0 0]';
ContactSpring3 = CalculateForce(Phi3,c);
Spring3 = [0 0 0 0 0 0 0 0 0 ContactSpring3 0 0]';


% Fextcall=Fextc-Fec;

Fextcall=Fextc-Fec+Spring1+Spring2+Spring3;                 % external forces minus elastic forces
Fextcsall=sparse(Fextcall);         % sparsity

y2sdot=MInvSparse*Fextcsall;        % accelerations (linear constraints are eliminated)
y2dot=full(y2sdot);                 % reconstruct full matrix from sparse matrix

% The first order differential equations (ODE) 
yddot(1:ndof)=y2dot(1:ndof);

