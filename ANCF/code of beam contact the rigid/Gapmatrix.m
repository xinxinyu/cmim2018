function PHI = Gapmatrix(t,u,v,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse,dt,mu,DIM,nc,Nk,ee,xi)

% build gap function PHI
% xi= linspace (-1,1,nc);
Ss1 = shapefunc_2322(H,-1,-1,DIM);
Phi1 = Ss1*ee(xloc(1,:));
PHI1 = Phi1(2);

for i = 1:nl
        xlock = xloc(i,:);
        eek=ee(xlock);  
        
for j = 2:nc
    XI = xi(j);
    ETA = -1;
    Ss2 = shapefunc_2322(H,XI,ETA,DIM); 
    Phi2 = Ss2*eek;
    m = 2*(j-1)-1: 2*(j-1);
    PHi2(m,1) = Phi2;
end
    n = 2*(nc-1)*i-2*nc+3:2*(nc-1)*i;
    PHIi2(n,1) = PHi2;
    k= (nc-1)*i-nc+2:(nc-1)*i;
    PHI2(k,1) = PHIi2(2*k);
end
PHI=[PHI1;
     PHI2]; 