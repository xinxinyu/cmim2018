% calculates right side of the ode
function ydot = eom_Dynam(t,y,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse)
% nx   DOFs 
% ndof DOFs after linear constraints

ee=zeros(nx,1);         % positions, systeemin kaikki koordinaatit, sek� lineaarisesti rajoitetut
eedot=zeros(nx,1);      % velocities, -"-

% positions
ee(bc)=y(1:ndof);
% velocities
eedot(bc)=y(ndof+1:2*ndof);  

if bcInd~=0
    ee(bcInd)=ee0bc;
    eedot(bcInd)=ee0dotbc;
end

ee=ee(:);
eedot=eedot(:);

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

Fextcall=Fextc-Fec;                 % external forces minus elastic forces
Fextcsall=sparse(Fextcall);         % sparsity

y1dot=eedot(bc);                    % velocities (linear constraints are eliminated)
y2sdot=MInvSparse*Fextcsall;        % accelerations (linear constraints are eliminated)
y2dot=full(y2sdot);                 % reconstruct full matrix from sparse matrix

% The first order differential equations (ODE) 
ydot(1:ndof)=y1dot(1:ndof); 
ydot(ndof+1:2*ndof)=y2dot(1:ndof);
ydot=ydot(:);
