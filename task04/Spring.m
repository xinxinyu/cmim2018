function dy = Spring(t,y,m,k)


dy=zeros(2,1);

dy(1)=y(2);
dy(2)=(-k/m)*y(1);

end