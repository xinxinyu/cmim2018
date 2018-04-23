function M = Mg(xloc,nl,nx,Element,rho,H,W,L)      
% Creates mass matrix of system using function MassM to get mass matrix of
% single element.
% For the students written by MKM

M=zeros(nx,nx);

% loop over all elements
for k = 1:nl
    xlock = xloc(k,:);
    %eek=ee(xlock);        
    Mk=MassM(Element,rho,H,W,L);
    M(xlock,xlock) = M(xlock,xlock) + Mk; 
end

