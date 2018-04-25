function [t,y] = RK4exp(tspan,y0,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MgcInvs)
% Explicit Runge-Kutta 4

% constant step size is taken from tspan.
deltat=tspan(2);
t_end=tspan(end);

step=0;
yk=y0;
for time=0:deltat:t_end
    step=step+1;  
   
    k1=deltat*eom_Dynam(time,yk,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MgcInvs);
    k2=deltat*eom_Dynam(time,yk+k1/2,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MgcInvs);
    k3=deltat*eom_Dynam(time,yk+k2/2,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MgcInvs);
    k4=deltat*eom_Dynam(time,yk+k3,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MgcInvs);    
    yk=yk+(1/6)*(k1+2*k2+2*k3+k4);

    t(step,1)=time;
    y(step,:)=yk'; 
end