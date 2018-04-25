function S = shapefunc_2322(H,xi,eta,DIM)

%     xi=XI(1);
%     eta=XI(2);


    Svec=[(xi + 1)^2/2 - (3*xi)/2 - 1/2;
        (H*eta)/2 + (H*eta*(xi + 1)^2)/4 - (3*H*eta*(xi + 1))/4;
        2*xi - (xi + 1)^2 + 2;
        H*eta*(xi + 1) - (H*eta*(xi + 1)^2)/2;
        (xi + 1)^2/2 - xi/2 - 1/2;
        (H*eta*(xi + 1)^2)/4 - (H*eta*(xi + 1))/4];
    
    % construct the matrix
    for ii=1:DIM
      for jj=1:6
          jj2=(ii-1)+(jj-1)*2+1;
          S(ii,jj2)=Svec(jj);         
      end
    end
end
