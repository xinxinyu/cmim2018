 % Code basis:
% For the course : Koneen simuloinnin eoj 2011 coded by MKM 2011
% Recoded and heavily commented for Xinxin's research by MKM 14.3.2018 

clc;
clear all;
close all;
format long;

% For integrator
%IntegFormu=24;
IntegFormu=300;                      % Choose integrator you wanna use
t_end=2.5;                             % end time of simulation
t0=0;                                % w.l.o.g. t0=0
deltat=0.0001;                       % only relevant for the output in matlab's integrators
tspan=(0:deltat:t_end)';

Element=2322;
n=16;                                 % number of elements
nc = 3;     %% number of contact points at each element
Nk = 1+n*(nc-1);  %% number of contacts points of the beam

% XI for shape function
xi= linspace (-1,1,nc);

% Tolerances - for Matlab integrators
RelTol=1e-3; %-3
AbsTol=1e-6; %-6

% Element definitions
% Define element properties

Element=2322;
ElemType='2D three-noded 12 dofs ANCF beam';  
ElemNodes=3;
DofsAtNode=4;
ElemDofs=ElemNodes*DofsAtNode;
DIM=2;                                    


% Geometry parameters etc.
% Element assembly
L=1;
W=1;              % beam width (z-direction)
H=0.1;              % beam height (y-direction)
rho=7850;           % density
% rho=2700;           % density
g=9.81;
% E=2.1E7;
E=2.1E11*0.001;
% E=70E9;
nu=0.3;
ks=1;
G=E/(2*(1+nu));
Gks=ks*G;
lambda=E*nu/((1+nu)*(1-2*nu));
mu = 0.3;% Friction coefficient    % Xinxin rectified this part

%%[P0,nloc] = linemesh_2322(n,L);
[P0,nloc] = linemesh_2322_Xinxinproblem(n,L);                                                                                                                                                                                                                                                                                                                                                                                                                                            

xloc=xlocAllANCF_2322(nloc);

% get number of nodes (nn), number of elements (nl) and total number of DOFs (nx)
[nn,~] = size(P0);      
[nl,~] = size(nloc);
nx = DofsAtNode*nn;  % number of degrees of freedom of elements

% node coordinates are formed over the entire structure P0:sta
for jj=1:nn,
    ee0((jj-1)*DofsAtNode+1:(jj-1)*DofsAtNode+DofsAtNode)=P0(jj,:);
end

ee0=ee0(:);
ee0dot=zeros(nx,1);  % The initial configuration, all zeros in our example,

% boundary conditions initialization
bc = true(1,nx);

%bcInd=xlocANCF_2322(1,(1:2));       % simply supprted at another end
bcInd=0;                             % no linear constraints

if bcInd~=0
    bc(bcInd)=0;                       % number of degrees of freedom of system after linear constraints
end

ndof = sum(bc); 
    
% initialiazation y0 that includes position and velocites after bc
% elimination
y0=zeros(2*ndof,1);             % y contains the positions and velocities
y0(1:ndof)=ee0(bc);             % initial positions
y0(ndof+1:2*ndof)=ee0dot(bc);   % initial velocities

if bcInd~=0
    nbc=length(bcInd);              % number of linear constraints  
    ee0bc=ee0(bcInd);                 % constrained dofs
    ee0dotbc=ee0dot(bcInd);
else
    nbc=0;              % number of linear constraints  
    ee0bc=0;                 % constrained dofs
    ee0dotbc=0;
end    

% OPTIONS FOR MATLAB INTEGRATORS:
options = odeset('AbsTol',AbsTol,'RelTol',RelTol);

Mg=Mg(xloc,nl,nx,Element,rho,H,W,L/n);      % Elements parameters as input
Mgc=Mg(bc,bc);
MgcInvs=sparse(inv(Mgc));

tic;
display(['Beginning with time integration, t = 0, t_end = ',num2str(t_end),' s'])
tdisp=0;
if IntegFormu == 45 % ode45
    display('Integrator: ode45 (Matlab), options:')
    options
    [t,y]=ode45(@eom_Dynam,tspan,y0,options,ElemDofs,rho,g,lambda,G,L/n,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MgcInvs);
elseif IntegFormu == 15
    display('Integrator: ode15s (Matlab), options:')
    options
    [t,y]=ode15s(@eom_Dynam,tspan,y0,options,ElemDofs,rho,g,lambda,G,L/n,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MgcInvs);
elseif IntegFormu == 23
    display('Integrator: ode23s (Matlab), options:')
    options
    [t,y]=ode23s(@eom_Dynam,tspan,y0,options,ElemDofs,rho,g,lambda,G,L/n,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MgcInvs);
elseif IntegFormu == 24
    display('Integrator: ode23t (Matlab), options:')
    options
    [t,y]=ode23t(@eom_Dynam,tspan,y0,options,ElemDofs,rho,g,lambda,G,L/n,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MgcInvs);
    
elseif IntegFormu == 100
    % % % % % Coded explicit RK4 ***************************************
    [t,y]=RK4exp(tspan,y0,ElemDofs,rho,g,lambda,G,L/n,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MgcInvs);    

elseif IntegFormu == 200
    % % % % % Coded Semi implicit Euler ***************************************
    u0 = y0(1:ndof);                             %% Iniatial coordinates
    v0 = y0(ndof+1:2*ndof);                      %% Initial velocities
    [t,u,v]=odeSemiFE(tspan, u0, v0,ElemDofs,rho,g,lambda,G,L/n,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MgcInvs);
    y=[u',v'];   
 elseif IntegFormu == 300 
    % % % % % Coded generalized version with generalized contact points***************************************
    u0 = y0(1:ndof);                             %% Iniatial coordinates
    v0 = y0(ndof+1:2*ndof);                      %% Initial velocities
    [t,u,v]=odeSemiFE_contact_generalize_contactpoint(tspan, u0, v0,ElemDofs,rho,g,lambda,G,L/n,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MgcInvs,mu,DIM,nc,Nk,xi);
    y=[u',v'];
 elseif IntegFormu == 400 
    % % % % % Coded generalized version with 5 contact points at each element  ***************************************
    u0 = y0(1:ndof);                             %% Iniatial coordinates
    v0 = y0(ndof+1:2*ndof);                      %% Initial velocities
    [t,u,v]=odeSemiFE_contact_generalize(tspan, u0, v0,ElemDofs,rho,g,lambda,G,L/n,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MgcInvs,mu,DIM);
    y=[u',v'];
     
 elseif IntegFormu == 500
    % % % % % Coded Semi implicit Euler ***************************************
    u0 = y0(1:ndof);                             %% Iniatial coordinates
    v0 = y0(ndof+1:2*ndof);                      %% Initial velocities
    [t,u,v]=odeSemiFE_contact(tspan, u0, v0,ElemDofs,rho,g,lambda,G,L/n,H,W,ee0bc,ee0dotbc,bcInd,xloc,nl,nx,ndof,bc,MgcInvs,mu);
    y=[u',v'];

else
    disp('Please select a valid integrator.  Error in file ''run_Dynam''.');
end

save Xinxindata16.mat u

CPUtime=toc;
disp(['CPU-time: ' num2str(CPUtime)]);


% % POST PROCESSING
% 
% Recover solution for full nodal coordinate vector
eedat=zeros(length(t),nx);

if bcInd~=0,
    for ii=1:length(t),  
        eedat(ii,bcInd)=ee0(bcInd)';       
    end
end
    
eedat(:,bc)=y(:,1:ndof);

% Recover solution for time derivative of the nodal coordinate vector
eedotdat=zeros(length(t),nx);
eedotdat(:,bc)=y(:,ndof+1:2*ndof);


% figure(1) % trajectory
set(gca, 'FontSize', 18, 'FontName','Times New Roman');
set(text, 'FontSize', 18, 'FontName','Times New Roman');
set(legend, 'FontSize', 18, 'FontName','Times New Roman');

res=eedat(:,(xlocANCF_2322(nn,[1:2])));

plot(res(:,1),res(:,2),'r-');
title('Trajectory of the free end of the beam')
xlabel('r_1')
ylabel('r_2')

% 
figure(2) % Vertical position of the free end over time
set(gca, 'FontSize', 18, 'FontName','Times New Roman');
set(text, 'FontSize', 18, 'FontName','Times New Roman');
set(legend, 'FontSize', 18, 'FontName','Times New Roman');

plot(t,res(:,2),'r-');
title('Vertical position of the free end over time')
xlabel('t')
ylabel('r_2')
display('Calculation completed.')
% display('Please press any button to continue with visualization.')
% pause


% %% Dynamic video 2d plot
% set(gcf,'color','w');
% filename='vid.avi';
% v = VideoWriter(filename);
% v.FrameRate=0.1/2/0.0001;  
% open(v);
%   
% for ii=1:10:size(eedat,1) %1:5000:size(eedat,1)
%     clf
% %     tic
% %     figure(3)
%     
%     title(['t = ',num2str(tspan(ii),'%5.2f')])
%     xlabel('x')  
%     ylabel('y')
%     
%     %axis equal
%     axis equal
%     xlim([-1.5,1.5])
%     ylim([-1,6]) 
%     grid on
%     hold on
%     
%         
%     for k = 1:nl  
%         nodes=nloc(k,1:3);
%       %  ee0k=eedat(1,xlocANCF_2322(nodes,[1:DofsAtNode]));
%         eek=eedat(ii,xlocANCF_2322(nodes,[1:DofsAtNode]));
%         run_plot_v005_beams(Element,H,W,L/n,E,G,lambda,1,eek,DIM,[1,1,1]);%rand(1,3)
%         %run_plot_v005_beams(Element,H,W,L/n,E,G,lambda,1,eek,DIM,[1,1,1]);
%         %run_plot_wire_2322(ee0k,eek,H,W,L,E,G,ks,Element,DIM)
%         
%     end    
%     %grid minor
% %     drawnow update expose
% %     timeElapsed = toc;
% %     pause(deltat-timeElapsed);
% hold off
%     writeVideo(v,getframe(gcf));
% end
% close(v);
% 
% %% figure plot
% ii=25000;
% 
%     set(gcf, 'PaperUnits', 'centimeters');
%     set(gcf, 'PaperPosition', [0 0 10 8]);
%     xlabel('Position X (m)','FontUnits','points','interpreter','latex','FontSize',10,'FontName','Times')
%     ylabel('Position Y (m)','FontUnits','points','interpreter','latex','FontSize',10,'FontName','Times')
%     title(['t = ',num2str(tspan(ii),'%5.2f')])
%     axis equal
%     xlim([-0.3,1])
%     ylim([-0.1,0.5]) 
%     grid on
%     hold on
%     
%     for k = 1:nl  
%         nodes=nloc(k,1:3);
%       %  ee0k=eedat(1,xlocANCF_2322(nodes,[1:DofsAtNode]));
%         eek=eedat(ii,xlocANCF_2322(nodes,[1:DofsAtNode]));
%         run_plot_v005_beams(Element,H,W,L/n,E,G,lambda,1,eek,DIM,[1,1,1]);%rand(1,3)
%         %run_plot_v005_beams(Element,H,W,L/n,E,G,lambda,1,eek,DIM,[1,1,1]);
%         %run_plot_wire_2322(ee0k,eek,H,W,L,E,G,ks,Element,DIM)
%         
%     end    
%     hold on
%     line([-0.3,1],[0,0],'color','black','LineWidth', 1)
%     
% hold off
%   print -depsc2 frame25000.eps









% plot(tspan',u(2,:),'k','LineWidth', 1.5)
% grid on
% xlabel('Time $t$ (s)','FontUnits','points','interpreter','latex','FontSize',12,'FontName','Times')
% ylabel('Displacement of the y direction of the first node (m) ','FontUnits','points','interpreter','latex','FontSize',12,'FontName','Times')
% %title('Displacement of the y direction of the first node','FontUnits','points','FontWeight','normal','FontSize',16,'FontName','Times')
% grid on
% print -depsc2 Displace_ydirection1.eps
% print -dpdf Displace_ydirection1.pdf
% hold off






