function [K,ff,Fext]=Kt(ee,Case,Mz,Fy,xloc,nx,nl,nn,Element,CrossSec,ElemDofs,DofsAtNode,E,G,nu,ks3,ks2,ks1,lambda,A,Iz,Iy,J,H,W,L);

% Computes tangential matrix for entire system using finite difference method
% also returns vectors of resultant forces and external forces

eps1=2*10^(-6);                 
hvec=L*eps1*ones(nx,1);

K=zeros(nx,nx);

% Loop over all dofs
for j=1:nx
        I_vec=zeros(nx,1);
        I_vec(j)=1;
        
        h=hvec(j); 
                               
        [f0,~]=fun(ee-h*I_vec,Case,Mz,Fy,xloc,nx,nl,nn,Element,CrossSec,ElemDofs,DofsAtNode,E,G,nu,ks1,ks2,ks3,lambda,A,Iz,Iy,J,H,W,L);
        [f1,~]=fun(ee+h*I_vec,Case,Mz,Fy,xloc,nx,nl,nn,Element,CrossSec,ElemDofs,DofsAtNode,E,G,nu,ks1,ks2,ks3,lambda,A,Iz,Iy,J,H,W,L);        
        
        K(:,j)=((f1)-(f0))/(2*h);
end

[ff,Fext] = fun(ee,Case,Mz,Fy,xloc,nx,nl,nn,Element,CrossSec,ElemDofs,DofsAtNode,E,G,nu,ks1,ks2,ks3,lambda,A,Iz,Iy,J,H,W,L);
