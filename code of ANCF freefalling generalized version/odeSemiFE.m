function [ T, u,v] = odeSemiFE(tspan, u0, v0,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse)

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
    v(:,i) = v(:,i-1) + dt*eom_Dynam_g(T(i-1), u(:,i-1),ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse)';
    u(:,i) = u(:,i-1) +dt*eom_Dynam_f(T(i-1), v(:,i),ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse)';  
end



% for i = 2:n
%     dt = T(i)-T(i-1);
%     v(:,i) = v(:,i-1) + dt*eom_Dynam_g(T(i-1), u(:,i-1),ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse);
%     u(:,i) = u(:,i-1) + dt*eom_Dynam_f(T(i-1), v(:,i),ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse);
% end





% time step
% dt = tspan(1);
% 
% T = 0:dt:tspan(2);
% 
% n = length(T);
% 
% % Y = zeros(length(y0), n);
% 
% Y(:, 1) = y0(:);
% 
% V(1:3, 1) = y0(4:6)';
% V(4:6, 1) = -1*(M^-1)*K*y0(1:3)';
% 
% for i = 2:n                                  
%     
%     V(:, i) = V(:, i-1) + dt * fun(T(i-1), V(:, i-1));
%     Y(:, i) = Y(:, i-1) + dt * V(:, i);
% end






