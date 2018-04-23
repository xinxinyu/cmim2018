function yddot = eom_Dynam_contact_g(t,u,v,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse,dt,mu)
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
    D = [1 0 -H/2 0 0 0 0 0 0 0 0 0 ;
         0 1 0 -H/2 0 0 0 0 0 0 0 0 ;
         0 0 0 0 1 0 -H/2 0 0 0 0 0 ;
         0 0 0 0 0 1 0 -H/2 0 0 0 0 ;
         0 0 0 0 0 0 0 0 1 0 -H/2 0 ;
         0 0 0 0 0 0 0 0 0 1 0 -H/2];                      %% D matrix for the contact    % Xinxin rectified this part
    D = D';
    N = D'*MInvSparse*D;                                   %% Matrix N for CCP formulation            %% Xinxin rectified it
    N=(N+N')/2;
    V = eedot;
    di0 = zeros(size(D,2),1);   %% di0 = [0 Phi1/dt 0 Phi2/dt 0 Phi3/dt]'  %% Xinxin rectified it

% shape function at the contact point P     %% Xinxin rectified it
    SP1 = [1 0 -H/2 0 0 0 0 0 0 0 0 0 ;
           0 1 0 -H/2 0 0 0 0 0 0 0 0 ]; 
    SP2 = [0 0 0 0 1 0 -H/2 0 0 0 0 0 ;
           0 0 0 0 0 1 0 -H/2 0 0 0 0 ];                    
    SP3 = [0 0 0 0 0 0 0 0 1 0 -H/2 0 ;
           0 0 0 0 0 0 0 0 0 1 0 -H/2 ];

    rP1 = SP1*u;                                              %% position value at the contact point P     %% Xinxin rectified it
    Phi1 = rP1(2);                                             %% gap function for CCP problem              %% Xinxin rectified it       
    rP2 = SP2*u; 
    Phi2 = rP2(2);
    rP3 = SP3*u;                                              %% position value at the contact point P     %% Xinxin rectified it
    Phi3 = rP3(2); 

% gap function for CCP problem              %% Xinxin rectified it   
    di0 = [0 Phi1/dt 0 Phi2/dt 0 Phi3/dt]';                                     %% di0 = [0 Phi1/dt 0 Phi2/dt 0 Phi3/dt]'    %% Xinxin rectified it
 
    di1 = di0+D'*(V+MInvSparse*Fextcf-MInvSparse*Fecf);         %% Xinxin rectified it
    p = di1;                                                    %% Xinxin rectified it

    A = kron(eye(3), [ 0 -1
                       1 -mu ] );    % contact/friction force inequalitiesb = [0 0]';  
    b = [0 0 0 0 0 0 ]';         % contact force is above to 0

%% output force    
    F = quadprog(N, p, A, b, [], [], [], [], [], optimoptions('quadprog',...
    'Algorithm','interior-point-convex','Display','off'));          % Eq (15)
    %myF(:,l) = F;
  if abs(F(2)) < 1e-6  
      F(1) = 0;             % friction force should be zero when normal force is zero
  end;
  if abs(F(4)) < 1e-6  
      F(3) = 0; 
  end;
  if abs(F(6)) < 1e-6  
      F(5) = 0; 
  end;

Ftoee=D*F;
Fextcsall=Fextcf-Fecf+Ftoee;                 % external forces minus elastic forces
% Fextcsall=sparse(Fextcall);         % sparsity

X=Ftoee/dt;
 

y2sdot=MInvSparse*Fextcsall;        % accelerations (linear constraints are eliminated)
y2dot=full(y2sdot);                 % reconstruct full matrix from sparse matrix

% The first order differential equations (ODE) 
yddot(1:ndof)=y2dot(1:ndof);
save Xinxindata.mat X