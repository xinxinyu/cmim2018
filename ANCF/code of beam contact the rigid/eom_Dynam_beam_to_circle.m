function yddot = eom_Dynam_beam_to_circle(t,u,v,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse,dt,mu,DIM,nc,Nk,xi,r)
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
% % construct D matrix
    D = Dmatrix_beam_to_circle(t,u,v,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse,dt,mu,DIM,nc,Nk,xi,ee);
    D = D.';
    N = D.'*MInvSparse*D;                                   %% Matrix N for CCP formulation            %% Xinxin rectified it
    N=(N+N')/2;
    V = eedot;

    
%% shape function at the contact point P     %% Xinxin rectified it
PHI = Gapmatrix_beam_to_element(t,u,v,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse,dt,mu,DIM,nc,Nk,ee,xi,r);

di0 = zeros(size(D,2),1);   %% di0 = [0 Phi1/dt 0 Phi2/dt 0 Phi3/dt]'  %% Xinxin rectified it 
di0 = zeros(2*Nk,1);
di0(2:2:end,1)= PHI/dt; %% di0 = [0 Phi1/dt 0 Phi2/dt 0 Phi3/dt]'    %% Xinxin rectified it
 

    di1 = di0+D'*(V+MInvSparse*Fextcf-MInvSparse*Fecf);         %% Xinxin rectified it
    p = di1;     %% Xinxin rectified it

    A = kron(eye(Nk), [ 0 -1
                       1 -mu ] );    % contact/friction force inequalitiesb = [0 0]';  
    b = zeros(2*Nk,1);         % contact force is above to 0

%% output force    
%       try
        F = quadprog(N, p, A, b, [], [], [], [], [], optimoptions('quadprog',...
            'Algorithm','interior-point-convex','Display','off'));          % Eq (15)
%       catch ex
%           disp('Error in quadprog');
%       end
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