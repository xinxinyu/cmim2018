function dv=SysPosition(t,v,m,k,c,A,omga)
t

dv=zeros(2,1); % velocity and accelation

dv(1)=v(2);
dv(2)=(-1/m)*(c*v(2)+k*v(1));
% dv(2)=(1/m)*(-k*v(1)-c*v(2)+sin(omga*t));

end