function F = FF_2322(in1,L,H,xi,eta)
%FF_2322
%    F = FF_2322(IN1,L,H,XI,ETA)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    06-Mar-2018 09:05:28

ee1 = in1(1,:);
ee2 = in1(2,:);
ee3 = in1(3,:);
ee4 = in1(4,:);
ee5 = in1(5,:);
ee6 = in1(6,:);
ee7 = in1(7,:);
ee8 = in1(8,:);
ee9 = in1(9,:);
ee10 = in1(10,:);
ee11 = in1(11,:);
ee12 = in1(12,:);
t2 = xi.*2.0;
t3 = t2+2.0;
t4 = xi+1.0;
t5 = t4.^2;
t6 = H.*t5.*(1.0./4.0);
t7 = 1.0./L;
t8 = H.*eta;
t9 = t8-H.*eta.*t3.*(1.0./2.0);
t10 = H.*eta.*(3.0./4.0);
t13 = H.*eta.*t3.*(1.0./4.0);
t11 = t10-t13;
t12 = H.*eta.*(1.0./4.0);
t14 = xi-1.0./2.0;
t15 = xi+1.0./2.0;
t16 = 1.0./H;
t17 = H.*(1.0./2.0);
t18 = t6+t17-H.*t4.*(3.0./4.0);
t19 = H.*t4;
t20 = t19-H.*t5.*(1.0./2.0);
F = reshape([t7.*(ee11.*(t12-H.*eta.*t3.*(1.0./4.0))+ee3.*t11-ee1.*t14-ee7.*t9-ee9.*t15+ee5.*xi.*2.0).*-2.0,t7.*(ee4.*t11-ee2.*t14-ee8.*t9-ee10.*t15+ee6.*xi.*2.0+ee12.*(t12-t13)).*-2.0,t16.*(ee3.*t18+ee7.*t20+ee11.*(t6-H.*t4.*(1.0./4.0))).*2.0,t16.*(ee4.*t18+ee8.*t20+ee12.*(t6-H.*t4.*(1.0./4.0))).*2.0],[2,2]);