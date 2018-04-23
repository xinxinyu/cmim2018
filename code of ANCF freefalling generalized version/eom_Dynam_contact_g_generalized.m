function yddot = eom_Dynam_contact_g_generalized(t,u,v,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse,dt,mu,DIM)
% nx   DOFs 
% ndof DOFs after linear constraints
%% CALCULATE the velocity of the element

ee=zeros(nx,1);         % positions, systeemin kaikki koordinaatit, sek� lineaarisesti rajoitetut
eedot=zeros(nx,1);      % velocities, -"-


% positions
ee(bc)=u(1:ndof);
eedot(bc)=v(1:ndof); 


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
% Calculation of the external and elastic impulses
Fextcf = Fextc*dt;                  %% external impulse     %% Xinxin rectified it
Fecf = Fec*dt;                      %% elastic impulse      %% Xinxin rectified it


%% contact force


xi1=-1;eta1=-1;
S1 = shapefunc_2322(H,xi1,eta1,DIM);
xi2=-0.5;eta2=-1;
S2 = shapefunc_2322(H,xi2,eta2,DIM);
xi3=0;eta3=-1;
S3 = shapefunc_2322(H,xi3,eta3,DIM);
xi4=0.5;eta4=-1;
S4 = shapefunc_2322(H,xi4,eta4,DIM);
xi5=1;eta5=-1;
S5 = shapefunc_2322(H,xi5,eta5,DIM);

Nk = 5+(nl-1)*4;   %mumber of contact points

D = zeros(2*Nk,ndof);
D(1:2,1:12) = S1;
for i = 1:nl   
xlock = xloc(i,:);
m = 8*i-7:8*i;  
D2(m,xlock) = [S2;
               S3;
               S4;
               S5];  
end
D(3:2*Nk,1:ndof)=D2;
    D = D.';
    N = D.'*MInvSparse*D;                                   %% Matrix N for CCP formulation            %% Xinxin rectified it
    N=(N+N')/2;
    V = eedot;
    di0 = zeros(size(D,2),1);   %% di0 = [0 Phi1/dt 0 Phi2/dt 0 Phi3/dt]'  %% Xinxin rectified it
    %%  D matrix is done

   
    
% shape function at the contact point P     %% Xinxin rectified it
Phi1 = S1*ee(xloc(1,:));
PHI1 = Phi1(2);
for k = 1:nl
        xlock = xloc(k,:);
        eek=ee(xlock);  
        
m = 8*k-7:8*k;  
Phi2(m,1) = [S2*eek;
             S3*eek;
             S4*eek;
             S5*eek];
n = 4*k-3:4*k;
PHI2(n,1) = Phi2(2*n);
end

PHI=[PHI1;
     PHI2]; 
di0 = zeros(2*Nk,1);
di0(2:2:end,1)= PHI/dt; %% di0 = [0 Phi1/dt 0 Phi2/dt 0 Phi3/dt]'    %% Xinxin rectified it
 

    di1 = di0+D'*(V+MInvSparse*Fextcf-MInvSparse*Fecf);         %% Xinxin rectified it
    p = di1;     %% Xinxin rectified it

    A = kron(eye(Nk), [ 0 -1
                       1 -mu ] );    % contact/friction force inequalitiesb = [0 0]';  
    b = zeros(2*Nk,1);         % contact force is above to 0

%% output force    
    F = quadprog(N, p, A, b, [], [], [], [], [], optimoptions('quadprog',...
    'Algorithm','interior-point-convex','Display','off'));          % Eq (15)
    %myF(:,l) = F;
     for i = 1 : Nk
        if abs(F(2*i)) < 1e-4  
            F(2*i-1) = 0; % friction force must be zero if normal force is zero
        end
    end
 
Ftoee=D*F;
Fextcsall=Fextcf-Fecf+Ftoee;                 % external forces minus elastic forces
% Fextcsall=sparse(Fextcall);         % sparsity

% X=Ftoee/dt;
 

y2sdot=MInvSparse*Fextcsall;        % accelerations (linear constraints are eliminated)
y2dot=full(y2sdot);                 % reconstruct full matrix from sparse matrix

% The first order differential equations (ODE) 
yddot(1:ndof)=y2dot(1:ndof);
% save Xinxindata.mat X