function r=kinem_v004_Beams(Element,H,W,L,E,G,lambda,ks,XI,ee,DIM)
% Fiilussa kuvataan elementin kinematiikka. Pyrit��n esitt�m��n muotofunktiot 
% xi-koordinaatistossa [-1..1]...ainakin k�ytett�ess� isoparametrisia
% ancf-elementtej�
% Fiilu on jatkoa shapefunc-koodeihin
% Recoded by MKM 2007

ee=ee(:);

if Element==3333   %  dof beam
    xi=XI(1);
    eta=XI(2);
    zeta=XI(3);    

    Svec=[1/2*xi*(-1+xi);
        1/4*H*xi*eta*(-1+xi);
        1/4*W*zeta*xi*(-1+xi);
        -(-1+xi)*(xi+1);
        -1/2*H*eta*(-1+xi)*(xi+1);
        -1/2*W*zeta*(-1+xi)*(xi+1);
        1/2*xi*(xi+1);
        1/4*H*xi*eta*(xi+1);
        1/4*W*zeta*xi*(xi+1)];    

    % Muotofunktiovektori matriisiksi
     for ii=1:3,
        for jj=1:9,
         jjj=((jj-1)*3+1)+(ii-1);
         S(ii,jjj)=Svec(jj);
        end
     end         
     
     r=S*ee;
     
elseif Element==2322   %  dof beam
    xi=XI(1);
    eta=XI(2);
    S = shapefunc_2322(H,xi,eta,DIM);
        
    r=S*ee;        
     
     
else
    disp('****** Elementill� ei indeksi� !! (virheilmoitus tiedostosta kinem_v004_Beam.m) ******');
end    

if DIM==2
        r(3)=XI(3);   % Yleisyyden takia, laitetaan Z-koordinaatiksi vaikkapa XI(3) 
end
