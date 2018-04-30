function M=MassM(Element,rho,H,W,L)


detF0=1/4*L*H;
  
        
M=rho*W*detF0*[  8/15,     0,          0,          0,  4/15,     0,          0,          0, -2/15,     0,          0,          0;
         0,  8/15,          0,          0,     0,  4/15,          0,          0,     0, -2/15,          0,          0;
         0,     0, (2*H^2)/45,          0,     0,     0,     H^2/45,          0,     0,     0,    -H^2/90,          0;
         0,     0,          0, (2*H^2)/45,     0,     0,          0,     H^2/45,     0,     0,          0,    -H^2/90;
      4/15,     0,          0,          0, 32/15,     0,          0,          0,  4/15,     0,          0,          0;
         0,  4/15,          0,          0,     0, 32/15,          0,          0,     0,  4/15,          0,          0;
     0,     0,     H^2/45,          0,     0,     0, (8*H^2)/45,          0,     0,     0,     H^2/45,          0;
     0,     0,          0,     H^2/45,     0,     0,          0, (8*H^2)/45,     0,     0,          0,     H^2/45;
 -2/15,     0,          0,          0,  4/15,     0,          0,          0,  8/15,     0,          0,          0;
     0, -2/15,          0,          0,     0,  4/15,          0,          0,     0,  8/15,          0,          0;
     0,     0,    -H^2/90,          0,     0,     0,     H^2/45,          0,     0,     0, (2*H^2)/45,          0;
     0,     0,          0,    -H^2/90,     0,     0,          0,     H^2/45,     0,     0,          0, (2*H^2)/45];


