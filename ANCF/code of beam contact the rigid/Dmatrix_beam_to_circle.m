function D = Dmatrix_beam_to_circle(t,u,v,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse,dt,mu,DIM,nc,Nk,ee,xi,r)


xcir=L/2;
ycir=0;
Ss1 = shapefunc_2322(H,-1,-1,DIM);
rr1 = Ss1*ee(xloc(1,:));
Xcon1 = rr1(1);   %% x coordinate of contact point
Ycon1 = rr1(2);   %% y coordinate of contact point
nn1 = [Xcon1-xcir;
       Ycon1-ycir];
Nn1 = nn1/norm(nn1);
Tt1 = [Nn1(2);
       -Nn1(1)];
A1 = [Tt1 Nn1];
DD1 = Ss1'*A1;

for k = 1:nl
        xlock = xloc(k,:);
        eek=ee(xlock);  

    for j = 2:nc
        XI = xi(j);
        ETA = -1;
        Ss2 = shapefunc_2322(H,XI,ETA,DIM); 
        rr2= Ss2*eek;
        Xcon2 = rr2(1);   %% x coordinate of contact point
        Ycon2 = rr2(2);   %% y coordinate of contact point
        nn2 = [Xcon2-xcir;
               Ycon2-ycir];
        Nn2 = nn2/norm(nn2);
        Tt2 = [Nn2(2);
              -Nn2(1)];
        A2 = [Tt2 Nn2];
        D2 = Ss2'*A2;

        m = 2*(j-1)-1: 2*(j-1);
        Dd2(1:12,m) = D2;
    end
end

for i = 1:nl
    xlock = xloc(i,:);
    n = 2*(nc-1)*i-2*nc+3:2*(nc-1)*i;
    DD2(xlock,n)=Dd2;
end

D = zeros(ndof,2*Nk);
D(1:12,1:2) = DD1;
D(1:ndof,3:2*Nk)=DD2;