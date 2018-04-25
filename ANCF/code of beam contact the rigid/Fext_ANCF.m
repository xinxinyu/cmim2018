function Fext=Fext_ANCF(Case,Element,DofsAtNode,xlocAll,nl,nn,nx,Mz,Fy,H,W,ee)

% Define vector of external forces

Fext = zeros(nx,1);

if Element==3333   %
    if Case==1
       Fext(xlocANCF_3333(nn,2)) = Fy;  % Beam tip load
    else
        disp('****** Check Case !! (Error in the file Fext_ANCF_v11.m) ******');
    end
elseif Element==2322   %
    if Case==1
       Fext(xlocANCF_2322(nn,2)) = Fy;  % Beam tip load
    else
        disp('****** Check Case !! (Error in the file Fext_ANCF_v11.m) ******');
    end    
    
else
    disp('****** Check element !! (Error in the file Fext_ANCF_v11.m) ******');
end