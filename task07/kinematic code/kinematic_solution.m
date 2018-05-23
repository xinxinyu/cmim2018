function [t,q,dq,ddq] = kinematic_solution(model_def,tspan,q0,dq0,ddq0,d1,d2)


q(:,1) = q0;
dq(:,1) = dq0;
ddq(:,1) = ddq0;
cnt = 0;
t=0;
tsim(1)=0; 
%% Solution values
%maximum number of iterations
maxiter=100;
%parameter to establish the convergenvce criteria
epsilon=0.0000001;

[ C, Cq, Cp, G ] = Kinematic_Constraint( q0,dq0, model_def,d1,d2,t);
% [ C, Cq, Cp, G ] = Constr4bar( q0,dq0,t );

for t = tspan
  
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
tmpdq = dq(:,cnt);
[ C, Cq, Cp, G ] = Kinematic_Constraint( tmp,tmpdq, model_def,d1,d2,t);
% [ C, Cq, Cp, G ] = Constr4bar( tmp,tmpdq,t );
tsim(cnt+1)=t;
end

%% Estimation of positions for next step
q(:,cnt+1)=q(:,cnt);
%Presentation and store of the result values
q;
nconverg;
% Velocity  
dq(:,cnt) = Cp;
dq(:,cnt+1) = dq(:,cnt);
% Accelerations
ddq(:,cnt) = G;
ddq(:,cnt+1)=dq(:,cnt);
end
end