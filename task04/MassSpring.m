function [dy,dv] = MassSpring(t,y,v,m,k)


% dy=0;
% dv=0;

dy=v;             %% velocity of the system
dv=(-k/m)*y;      %% accelation of the system

end