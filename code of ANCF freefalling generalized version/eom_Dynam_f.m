% calculates right side of the ode
function ydot = eom_Dynam_f(t,v,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse,D)
% nx   DOFs 
%% CALCULATE the position of the element
% ndof DOFs after linear constraints

 eedot=zeros(nx,1);      % velocities, -"-

%velocities
eedot(bc)=v(1:ndof); 

if bcInd~=0
    eedot(bcInd)=ee0dotbc;
end

eedot=eedot(:);

y1dot=eedot(bc);                    % velocities (linear constraints are eliminated)


% The first order differential equations (ODE) 
ydot(1:ndof)=y1dot(1:ndof); 
