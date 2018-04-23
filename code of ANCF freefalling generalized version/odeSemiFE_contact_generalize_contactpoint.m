function [ T, u,v] = odeSemiFE_contact_generalize_contactpoint(tspan, u0, v0,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse,mu,DIM,nc,Nk,xi)

%   Simple integration using Semi-implicit Euler Method
%   fun - function handle with interface fun(t, y)
%   tspan - two element vector with dt and tend
%   y0 - initial conditions

T = tspan;
n = length(T);
u = zeros(length(u0),n);
v = zeros(length(v0),n);
u(:,1) = u0;
v(:,1) = v0;

for i = 2:n
    dt = T(i)-T(i-1);
    v(:,i) = v(:,i-1) + eom_Dynam_contact_g_generalized_contactpoint(T(i-1), u(:,i-1),v(:,i-1),ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse,dt,mu,DIM,nc,Nk,xi)';
    u(:,i) = u(:,i-1) +dt*eom_Dynam_contact_f(T(i-1), v(:,i),ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse)';  
end