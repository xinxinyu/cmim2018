% Basic static analysis code for beam elements
% Based on code by Marko Matikainen
% Modified by Vesa-Ville Hurskainen
% for simplicity, implements only 1 type of element (3D ANCF beam)

% 05032018 by MKM: Constinuum mechanics based three node element (2322) for Xinxin's usage.
% The elastic forces of the 2322 derived in Matlab.

clc;
clear
close all
format long;

Ffigplot=1;     % Flag for plotting (1 = plot deformation, 0 = don't)
scale=1;        % Scaling factor for displacement plot.

% Choose element: only 3333 (Nachbagauer quadratic 3D ANCF beam)
%Element=3333;
Element=2322;
nmesh=[1,2,4,8,16];     % Vector of meshes (element numbers) to be computed

% Choose element cross section: only 1 = rectangular
CrossSec=1;

% Select problem for computation
Case=1;

% Define element properties
if Element==3333
    ElemType='3D three-noded 27 dofs ANCF beam';  
    ElemNodes=3;
    DofsAtNode=9;
    ElemDofs=ElemNodes*DofsAtNode;
    DIM=3;  
elseif Element==2322 % Xinxin: change it for 2322 
    ElemType='2D three-noded 12 dofs ANCF beam';  
    ElemNodes=3;
    DofsAtNode=4;
    ElemDofs=ElemNodes*DofsAtNode;
    DIM=2;                                    
    
    
else
    disp('****** Choose a valid element !! (statictest_ANCFbeam.m) ******');  
end

% Define problem
if Case==1 % Large deformation cantilever problem
    H=0.5;
    W=0.1;
    L=2;     
    E=2.07e11; 
    nu=0.3;
    ks3=10*(1+nu)/(12+11*nu);       % shear correction factor (around z-axis)
    ks1=1;                          % shear correction factor (around x-axis)
    ks2=ks3;                        % shear correction factor (around y-axis)

    Fy=-5e8*H^3;  % Beam tip load
    Mz=0;         % Moment load (not used atm!)
    
    G=E/(2*(1+nu));% shear modulus
    lambda=E*nu/((1+nu)*(1-2*nu)); % construct the matrix E, E is the matrix of elastic constants of the material
     
    % Cross-sectional parameters for structural elements (not used atm!)
    A=W*H;
    Iz=W*H^3/12;
    Iy=H*W^3/12;
    J = 0.291*H*W^3;
    
else
     disp(['Choose a valid case !'])
end

Case1 = [];

for n=nmesh         % Loop over all defined meshes: nmesh=[1 ....])
    clear P P0 xloc u0 uu ee bc f1 K ff ffcs Kc Kcs;    % clear previous mesh definitions
 
    if Element==3333
        [P0,nloc] = linemesh3nodes_v10_ANCF_3333(n,L);    
        xloc=xlocAllANCF_3333(nloc);                   
        
        % get number of nodes (nn), number of elements (nl) and total number of DOFs (nx)
        [nn,nodedof] = size(P0);
        [nl,m] = size(nloc);
        nx = 9*nn;          % dofs (no constraints) % Xinxin: change it for 2322 
        
        % draw system
        u0 = zeros(nx,1);
        uu = u0;
        
        % create global vector of nodal coordinates
        for jj=1:nn,
            ee((jj-1)*9+1:(jj-1)*9+9)=P0(jj,:); % Xinxin: change it for 2322 
        end
        
        % Define vector of linear constraints
        bc = logical(ones(1,nx));
        bc(xlocANCF_3333(1,[1:3,4,6,7,8]))=0;       % don't constrain cross-sectional deformation
        %bc(xlocANCF_3333(1,[1:9]))=0;              % constrain all DOFs
        
    elseif Element==2322  
        [P0,nloc] = linemesh_2322(n,L);  
        xloc=xlocAllANCF_2322(nloc);                    
        
        % get number of nodes (nn), number of elements (nl) and total number of DOFs (nx)
        [nn,nodedof] = size(P0);
        [nl,m] = size(nloc);
        nx = 4*nn;          % dofs (no constraints) % Xinxin: change it for 2322 
        
        % draw system
        u0 = zeros(nx,1);
        uu = u0;
        
        % create global vector of nodal coordinates
        for jj=1:nn,
            ee((jj-1)*4+1:(jj-1)*4+4)=P0(jj,:); % Xinxin: change it for 2322 
        end
            
        
        % Define vector of linear constraints
        bc = logical(ones(1,nx));
        bc(xlocANCF_2322(1,[1:4]))=0;       % don't constrain cross-sectional deformation, 
               
    
    else
        disp('****** Choose a valid element !! (statictest_ANCFbeam.m) ******');  
    end
    
    % Define initial position
    ee=ee(:);
    e0=ee;
    
    ndof = sum(bc);     % Number of unconstrained DOFs
            
    imax=50;            % Maximum number of iterations for Newton's method
    
    titertot=0;
    
    % START NEWTON'S METHOD
    for ii=1:imax,
        
        tic;
        [K,ff,Fext] = Kt(ee,Case,Mz,Fy,xloc,nx,nl,nn,Element,CrossSec,ElemDofs,DofsAtNode,E,G,nu,ks3,ks2,ks1,lambda,A,Iz,Iy,J,H,W,L/n);   
        % Eliminate linear constraints
        Kc = K(bc,bc);
        ffc=ff(bc); 
        
        deltaf=ffc/norm(Fext(bc));  % compute residual

        uuc = -Kc\ffc;          % compute displacements
        uu(bc) = uu(bc)+uuc;    % add displacement to initial condition
        ee(bc) = ee(bc)+uuc;
        
        titer=toc;
        titertot=titertot+titer;
        disp(['Iteration: ' num2str(ii) '  ,CPU-time: ' num2str(titer)]);
        
        Re=10^(-4);             % Stopping criterion for residual

        if abs(deltaf)<Re*ones(ndof,1)
            disp(['Solution is found by ' num2str(ii) ' iterations. Total CPU-time: ' num2str(titertot)])
            break

        elseif ii==imax 
            disp(['The solution is not found. The maximum number of iterations is reached. Total CPU-time: ' num2str(titertot)])
        else     
            disp(['wait...'])
        end              
 

    end
    % END NEWTON
      
    % Pick nodal displacements from result vector
    if Element==3333 % 
        ux = uu(xlocANCF_3333(nn,[1])); 
        uy = uu(xlocANCF_3333(nn,[2]));
        uz = uu(xlocANCF_3333(nn,[3]));
    elseif Element==2322 
        ux = uu(xlocANCF_2322(nn,[1])); 
        uy = uu(xlocANCF_2322(nn,[2]));        
        
    else
        disp('****** Choose a valid element !! (statictest_ANCFbeam.m) ******');  
    end  
    
    if DIM == 3,
        Case1 = [Case1; n ndof ux uy uz]; 
    else
        Case1 = [Case1; n ndof ux uy 0];
    end    
        
end


% POST PROCESSING

disp(sprintf('Nonlinear static test, Case  %g, Elem = %g, E = %g, nu = %g, ks=%g,H = %g, Fy=%g', Case, Element,E,nu,ks3,H,Fy))
disp('& Mesh & DOFs & ux & uy & uz \\)')
for k=1:size(Case1,1) 
  disp(sprintf('& %4d & %4d  & %10.6f & %10.6f & %10.6f\\\\',Case1(k,[1:5])))
end


% Visualization of displacement
if Ffigplot==1,
figure(1);
set(gca, 'FontSize', [20], 'FontName','Times New Roman');
set(text, 'FontSize', [20], 'FontName','Times New Roman'); 

xlabel('{\it{X}} [m]','FontName','Times New Roman','FontSize',[20]),ylabel('{\it{Y}} [m]','FontName','Times New Roman','FontSize',[20]),zlabel('Z [m]','FontName','Times New Roman','FontSize',[20]);

title(['Displacements: X:' ,num2str(ux) ,'  Y:', num2str(uy)],'FontName','Times New Roman','FontSize',[20]); 

for k = 1:nl
    if Element==3333
        nodes=nloc(k,1:3);
        e0k=e0(xlocANCF_3333(nodes,[1:DofsAtNode]));
        eek=ee(xlocANCF_3333(nodes,[1:DofsAtNode]));
    elseif Element==2322
        nodes=nloc(k,1:3);
        e0k=e0(xlocANCF_2322(nodes,[1:DofsAtNode]));
        eek=ee(xlocANCF_2322(nodes,[1:DofsAtNode]));    
        
        
    else
        disp('****** Choose valid element !! (statictest_ANCFbeam.m) ******');  
    end               
    run_plot_v005_beams(Element,H,W,L/n,E,G,lambda,ks3,eek,DIM,[1,1,1])
end
end








