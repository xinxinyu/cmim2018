function dv=SysPosition(t,v,p)

m=p.m;
k=p.k;
c=p.c;

dv=zeros(2,1);


dv(1)=v(2);
dv(2)=(-1/m)*(c*v(2)+k*v(1));


end