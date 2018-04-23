function [P,nloc] = linemesh_2322(n,L)

% make a linemesh with three node ANCF beam element 2322
% Coded by MKM
% ..v10 is improved version compared to the older linemesh code


% generate nodal coordinates in xy-plane
dx = L/n;
nodes=n*3-(n-1);

% Geospace for spacing
%xk=geospace(0,a,nodes,0.95)';

xk=geospace(0,L,nodes,1)';
yk=zeros(1,nodes)';

nullmat = zeros(nodes,2);
onesvec = ones(nodes,1); %% ones-Create array of all ones

%dr1dx=nullmat;
%dr1dx(:,1)=onesvec;

dr1dy=nullmat;
dr1dy(:,2)=onesvec;

%dr1dz=nullmat;
%dr1dz(:,3)=onesvec;


P = [];
Pk = [xk yk dr1dy];  
P = [P; Pk]; 

% generate elememnt connectivity
nloc = [];

for k = 1:n
    loc = [(k-1)*2+1 (k-1)*2+2 (k-1)*2+3];
    nloc = [nloc; loc];
end
