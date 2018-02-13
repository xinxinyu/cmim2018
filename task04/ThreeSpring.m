function dy = ThreeSpring(t,y,M,K)


dy=zeros(6,1);

dy(1:3)=y(4:6);

dy(4:6)=-(M^-1)*K*y(1:3);

end