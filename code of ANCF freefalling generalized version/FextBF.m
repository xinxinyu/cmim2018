function Fext=FextBF(rho,g,L,H,W)
% contribution to external forces by body forces
% For element 2322

if g~=0
   Smint=[2/3,0, 0, 0,8/3,0, 0, 0, 2/3,0, 0, 0;
          0, 2/3, 0, 0,   0, 8/3, 0, 0,   0, 2/3, 0, 0];
    b=[0; -rho*g];              % body forces due to gravity       
    detF0=L*H/4;
    Fext=b'*Smint*W*detF0;    
else
    Fext=zeros(2*4,1);
end

Fext=Fext(:);
 