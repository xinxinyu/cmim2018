function run_plot_v005_beams(Element,H,W,L,E,G,lambda,ks,ee,DIM,color)
% Plottausrutiini elementi kinematiikalle. Pelk�st��n B12

%L�ht� ja loppupisteet elementtikoordinaatistossa
x1=-1;
x2=1;

% y1=-H/2;
% y2=H/2;

y1=-1;
y2=1;

z1=-1;
z2=1;

x_scale=1/4;
y_scale=1/2;
z_scale=1/2;

% x_scale=1;
% y_scale=1;
% z_scale=1;

%Elementin mitat ja gridin jako x,y ja z suunnassa
Elem_dim=[x1,x2,x_scale;
          y1,y2,y_scale;
          z1,z2,z_scale];
      
scale_A=Elem_dim(1,3);
a=Elem_dim(1,1);
b=Elem_dim(1,2);
scale_B=Elem_dim(2,3);
c=Elem_dim(2,1);
d=Elem_dim(2,2);

e=Elem_dim(3,2);
f=Elem_dim(3,2);

%Suoraklulm s�rm tahkot
cubic=[Elem_dim(1,1), Elem_dim(1,2), Elem_dim(1,3), Elem_dim(2,1), Elem_dim(2,2), Elem_dim(2,3), Elem_dim(3,1);
       Elem_dim(1,1), Elem_dim(1,2), Elem_dim(1,3), Elem_dim(2,1), Elem_dim(2,2), Elem_dim(2,3), Elem_dim(3,2);
       Elem_dim(3,1), Elem_dim(3,2), Elem_dim(3,3), Elem_dim(2,1), Elem_dim(2,2), Elem_dim(2,3), Elem_dim(1,1);
       Elem_dim(3,1), Elem_dim(3,2), Elem_dim(3,3), Elem_dim(2,1), Elem_dim(2,2), Elem_dim(2,3), Elem_dim(1,2);
       Elem_dim(1,1), Elem_dim(1,2), Elem_dim(1,3), Elem_dim(3,1), Elem_dim(3,2), Elem_dim(3,3), Elem_dim(2,1);
       Elem_dim(1,1), Elem_dim(1,2), Elem_dim(1,3), Elem_dim(3,1), Elem_dim(3,2), Elem_dim(3,3), Elem_dim(2,2)];

%figure(1)
%view(20,12)
axis equal
%axis([-1.5 2.2 -1.2 1.2 -1 1]);
xlabel('X','Fontsize',12);
ylabel('Y','Fontsize',12);
zlabel('Z','Fontsize',12);
   
for cc=1:6
       a=cubic(cc,1);
       b=cubic(cc,2);
       scale_A=cubic(cc,3);
       c=cubic(cc,4);
       d=cubic(cc,5);
       scale_B=cubic(cc,6);
       e=cubic(cc,7);
       
       [AA,BB,CC,i_end,j_end]=func_face(a,b,c,d,e,scale_A,scale_B);
   
       for i=1:(i_end)
           for j=1:(j_end)

               patch_coord_A=[2, 2];
               patch_coord_B=[2, 2];
            
               patch_coord_A=AA(i:1+i,j:j+1);
               patch_coord_B=BB(i:1+i,j:j+1);
               patch_coord_C=CC(i:1+i,j:j+1);
       
               patch_coord_vec_A(1)=patch_coord_A(1,1);
               patch_coord_vec_A(2)=patch_coord_A(1,2);
               patch_coord_vec_A(3)=patch_coord_A(2,2);
               patch_coord_vec_A(4)=patch_coord_A(2,1);
       
               patch_coord_vec_B(1)=patch_coord_B(1,1);
               patch_coord_vec_B(2)=patch_coord_B(1,2);
               patch_coord_vec_B(3)=patch_coord_B(2,2);
               patch_coord_vec_B(4)=patch_coord_B(2,1);
       
               patch_coord_vec_C(1)=patch_coord_C(1,1);
               patch_coord_vec_C(2)=patch_coord_C(1,2);
               patch_coord_vec_C(3)=patch_coord_C(2,2);
               patch_coord_vec_C(4)=patch_coord_C(2,1);
       
 % Seuraava(kin) osio tehty typer�sti, koodaa joskus paremmin k�ytt�m�ll� aliohjelmia                            
               if cc==1 | cc == 2
                   for ii=1:4
                       XI(1)=patch_coord_vec_A(ii);
                       XI(2)=patch_coord_vec_B(ii);
                       XI(3)=patch_coord_vec_C(ii);
                       %r=kinem_v003(Element,H,W,L,E,G,ks,XI,ee,DIM);
                       r=kinem_v004_Beams(Element,H,W,L,E,G,lambda,ks,XI,ee,DIM); 
                       %[r]=shapefunc(H,B,L,x,y,z,ee);
                       AA_def(ii)=r(1);
                       BB_def(ii)=r(2);
                       CC_def(ii)=r(3);
                   end
                   patch(AA_def(1:4)',BB_def(1:4)',CC_def(1:4)', color);
                   % undef.
                   %patch(patch_coord_vec_A(1:4)', patch_coord_vec_B(1:4)', patch_coord_vec_C(1:4)', 0.5)
                   
               elseif cc == 3 | cc == 4
                   for ii=1:4
                       XI(3)=patch_coord_vec_A(ii);
                       XI(2)=patch_coord_vec_B(ii);
                       XI(1)=patch_coord_vec_C(ii);
                       %[r]=shapefunc(H,B,L,x,y,z,ee);
                       %r=kinem_v003(Element,H,W,L,E,G,ks,XI,ee,DIM);
                       r=kinem_v004_Beams(Element,H,W,L,E,G,lambda,ks,XI,ee,DIM);
                       AA_def(ii)=r(3);
                       BB_def(ii)=r(2);
                       CC_def(ii)=r(1);
                   end
                   
                   patch(CC_def(1:4)',BB_def(1:4)',AA_def(1:4)', color);
                   %undef.
                   %patch(patch_coord_vec_C(1:4)', patch_coord_vec_B(1:4)', patch_coord_vec_A(1:4)', 0.5)
                   
               else 
                   for ii=1:4
                       XI(1)=patch_coord_vec_A(ii);
                       XI(3)=patch_coord_vec_B(ii);
                       XI(2)=patch_coord_vec_C(ii);
                       %[r]=shapefunc(H,B,L,x,y,z,ee);
                       %r=kinem_v003(Element,H,W,L,E,G,ks,XI,ee,DIM);
                       r=kinem_v004_Beams(Element,H,W,L,E,G,lambda,ks,XI,ee,DIM);
                       AA_def(ii)=r(1);
                       BB_def(ii)=r(3);
                       CC_def(ii)=r(2);
                   end
                   patch(AA_def(1:4)',CC_def(1:4)',BB_def(1:4)', color);
                   %undef.
                   %patch(patch_coord_vec_A(1:4)', patch_coord_vec_C(1:4)', patch_coord_vec_B(1:4)', 0.5)
               end
               
            end
       end
end

view3d rot
