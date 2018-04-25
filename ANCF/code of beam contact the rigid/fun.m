function [ff,Fext] = fun(ee,Case,Mz,Fy,xloc,nx,nl,nn,Element,CrossSec,ElemDofs,DofsAtNode,E,G,nu,ks1,ks2,ks3,lambda,A,Iz,Iy,J,H,W,L)

% Forms force equation ff=Fe-Fext  (elastic forces - external forces)

Fe=zeros(nx,1);

% Loop over elements
for k = 1:nl
    xlock = xloc(k,:);
    eek=ee(xlock);
    
    if Element == 3333 % Nachbagauer 3D quadratic element
         if CrossSec==1,    
         % rectangular cross-section, enhanced continuum mechanics approach
         Fek=Integ_Fe_3333(ElemDofs,E,nu,G,ks2,ks3,H,W,L,eek)';  
            
         else
             sprintf('Error in function: fun!! Choose a valid cross section')   
         end
         
    elseif Element == 2322 % Nachbagauer 2D quadratic element
         if CrossSec==1,    
         % rectangular cross-section
         Fek=Fe_2322(ElemDofs,lambda,G,L,H,W,eek)'; 
         
         
            
         else
             sprintf('Error in function: fun!! Choose a valid cross section')   
         end     

    else
        sprintf('Error in function: fun!! Choose a valid element')   
    end

    Fe(xlock)=Fe(xlock)+Fek;

end


                   
Fext=Fext_ANCF(Case,Element,DofsAtNode,xloc,nl,nn,nx,Mz,Fy,H,W,ee);



ff=Fe-Fext;
ff=ff(:);