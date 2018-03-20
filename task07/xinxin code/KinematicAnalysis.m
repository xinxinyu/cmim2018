%% initial values of the system
q(:,1)=q0;
dq(:,1)=dq0;
ddq(:,1)=ddq0;
tsim(1)=tini; 

for t = tini:dt:tend
  
cnt = cnt+1; %Column counter to form the result matrix.
nloop=0;
nconverg=0;

while nloop<maxiter

%% Newton Difference
deltaq = -Cq^(-1)*C;
% Coordinate update
q(:,cnt)=q(:,cnt)+deltaq;

        nloop=nloop+1;
        nconverg=nconverg+1;    
 %%%% Converging criteria (error < epsilon)
maxdx=max(abs(deltaq));
        if maxdx<epsilon
            nloop=maxiter;
        end
 % constrains update       
tmp=q(:,cnt);
C = ConsEqua4Bar(tmp,t);
% % jacobian update
Cq = Jaco4Bar(tmp);
tsim(cnt+1)=t;
end

%% Estimation of positions for next step
q(:,cnt+1)=q(:,cnt);
%Presentation and store of the result values
q;
nconverg;
%Determination of the velocities (dq)
Ct = ConsTime(tmp,t);
% Velocity  
dq(:,cnt) = Cq^(-1)*-Ct;
dq(:,cnt+1) = dq(:,cnt);

%% Determination of the accelerations
tmpdq = dq(:,cnt);
Ctt = ConsTimeTime(tmpdq,t); 
Cqqq = Consqq(tmp,tmpdq,t);
Cqt = Consqt(tmpdq,t);
 
ddq(:,cnt)=Cq^(-1)*(-Ctt-Cqqq-Cqt);
ddq(:,cnt+1)=dq(:,cnt);
end


