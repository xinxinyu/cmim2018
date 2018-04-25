function [ T, u, v] = odeSemiFE_onepoint(tspan, u0, v0,ElemDofs,rho,g,lambda,G,L,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MInvSparse,mu)
% nx   DOFs 
% ndof DOFs after linear constraints
T = tspan;

n = length(T);
u = zeros(length(u0),n);
v = zeros(length(v0),n);
u(:,1) = u0;
v(:,1) = v0;
l = 1;
for l = 1:n-1
    dt = T(l+1)-T(l);                                      %% time step
    D = [0 0 0 0 1 0 -H/2 0 0 0 0 0 ;
         0 0 0 0 0 1 0 -H/2 0 0 0 0 ];                     %% D matrix for the contact    % Xinxin rectified this part
    D = D';
    N = D'*MInvSparse*D;                                   %% Matrix N for CCP formulation            %% Xinxin rectified it
    N=(N+N')/2;
    V = v(:,l);
    
    di0 = zeros(size(D,2),1);                              %% di0 = [0 Phi1/dt 0 Phi2/dt 0 Phi3/dt]'  %% Xinxin rectified it
    SP = [0 0 0 0 1 0 -H/2 0 0 0 0 0 ;
          0 0 0 0 0 1 0 -H/2 0 0 0 0 ];                    %% shape function at the contact point P     %% Xinxin rectified it
    rP(:,l) = SP*u(:,l);                                   %% position value at the contact point P     %% Xinxin rectified it
    Phi = rP(2,l);                                         %% gap function for CCP problem              %% Xinxin rectified it   
    di0 = [0 Phi/dt]';                                     %% di0 = [0 Phi1/dt 0 Phi2/dt 0 Phi3/dt]'    %% Xinxin rectified it
        
ee=zeros(nx,1);         % positions, systeemin kaikki koordinaatit, sek� lineaarisesti rajoitetut
eedot=zeros(nx,1);      % velocities, -"-

% positions
ee(bc)=u(1:ndof);
% velocities
eedot(bc)=v(1:ndof); 

if bcInd~=0
    ee(bcInd)=ee0bc;
    eedot(bcInd)=ee0dotbc;
end

ee=ee(:);
eedot=eedot(:);

%% Calculation about external and elastic forces
Fe=zeros(nx,1);
Fext=zeros(nx,1);

% loop over the elements - calculation of elastic and external forces
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
Fextc=Fext(bc);
Fec=Fe(bc);
%% Calculation of the external and elastic impulses
Fextcf = Fextc*dt;                  %% external impulse     %% Xinxin rectified it
Fecf = Fec*dt;                      %% elastic impulse      %% Xinxin rectified it

di1 = di0+D'*(V+MInvSparse*Fextcf+MInvSparse*Fecf);         %% Xinxin rectified it
p = di1;                                                    %% Xinxin rectified it

A = [ 0 -1
      1 -mu ];       % contact/friction force inequalities
b = [0 0]';  


%% output force    
    F(:,l) = quadprog(N, p, A, b, [], [], [], [], [], optimoptions('quadprog',...
    'Algorithm','interior-point-convex','Display','off'));          % Eq (15)
    %myF(:,l) = F;
  if abs(F(2)) < 1e-6  % friction force should be zero when normal force is zero
      F(1) = 0; 
  end
  
    v(:,l+1) = v(:,l) + (MInvSparse*Fextcf+MInvSparse*Fecf + MInvSparse*D*F(:,l));   % Eq 13
    u(:,l+1) = u(:,l) + v(:,l+1)*dt;          % Eq 9a

end

