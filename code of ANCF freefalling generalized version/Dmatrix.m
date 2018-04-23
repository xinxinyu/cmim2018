function D = Dmatrix(t,u,v,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse,dt,mu,DIM,nc,Nk,xi)

% construct D matrix
% xi= linspace (-1,1,nc);
Ss1 = shapefunc_2322(H,-1,-1,DIM);
D1 = Ss1;

for j = 2:nc
    XI = xi(j);
    ETA = -1;
    Ss2 = shapefunc_2322(H,XI,ETA,DIM); 
    m = 2*(j-1)-1: 2*(j-1);
    SS2(m,1:12) = Ss2;
end

for i = 1:nl
    xlock = xloc(i,:);
    n = 2*(nc-1)*i-2*nc+3:2*(nc-1)*i;
    D2(n,xlock)=SS2;
end

D = zeros(2*Nk,ndof);
D(1:2,1:12) = D1;
D(3:2*Nk,1:ndof)=D2;