function [AA,BB,CC,i_end,j_end]=func_face(a,b,c,d,e,scale_A,scale_B)
% Written by MKM LUT/IMVe 2004

if a<b & a*b ~= 0
    j_end=(b-a)/scale_A;
    ab=b;
elseif a>b & a*b ~= 0
    j_end=(a-b)/scale_A;
    scale_A=-scale_A
    ab=-a;
elseif (a == 0 | b ==0) & a<b
    j_end=(b-a)/scale_A;
    ab=-a;
elseif (a == 0 | b == 0) & a>b
    j_end=(a-b)/scale_A;
    scale_A=-scale_A;
    ab=0;
else
    disp(['VIRHE !!! PARAMETRIT YHTÄSUURIA !!!'])
    %break      
end
    
if c<d & c*d ~= 0
    i_end=(d-c)/scale_B;
    cd=d;
elseif c>d & c*d ~= 0
    i_end=(c-d)/scale_B;
    scale_B=-scale_B
    cd=-c;
elseif (c == 0 | d ==0) & c<d
    i_end=(d-c)/scale_B;
    cd=-c;
elseif (c == 0 | d == 0) & c>d
    i_end=(c-d)/scale_B;
    scale_B=-scale_B;
    cd=0;
else
    disp(['VIRHE !!! PARAMETRIT YHTÄSUURIA !!!'])
    %break    
end

AA=[i_end+1,j_end+1];

for i=1:(i_end+1)
   for j=1:(j_end+1)
       AA(i,j)=(j-1)*scale_A-ab;
   end
end

for i=1:(i_end+1)
   for j=1:(j_end+1)
       BB(i,j)=(i-1)*scale_B-cd;
   end
end

CC=ones(i_end+1,j_end+1)*e;

