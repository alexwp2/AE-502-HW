function [dydt,t,r] = GaussPlanetaryEqs (t,y,y_0,mu,J2,R_earth,x,period)
% r = x(round(t)+1);
% y_0 = [a_0, i_0, e_0, w_0, omega_0, M_0];
n = sqrt(mu/(y(1))^3);
t_0 = y_0(6)/n;
% t_new = rem(t,period);
M = rem(abs(y(6) + n*(t-t_0)),2*pi);
E = M + y(3)*sin(M) + y(3)^2/2*sin(2*M) + y(3)^3/8*(3*sin(3*M) - sin(M)) + y(3)^4/6*(2*sin(4*M) - sin(2*M));
v = 2*atan(sqrt((1+y(3)))*tan(E/2)/sqrt(1-y(3)));
% if v >= pi
%     v = v - pi;
% elseif v <= -pi
%     v = v + pi;
% end
p = y(1)*(1 - y(3)^2);
u = v + y(4);
r = p/(1+y(3)*cos(v)); 
P2 = 0.5*(3*(sin(u))^2*(sin(y(2)))^2 - 1);
dP2 = (-2*sin(u)*sin(y(2))*P2 + 2*sin(u)*sin(y(2)));
% Fr = 3*mu*J2*R_earth^2*(3*sin(u)^2*sin(y(2))^2 - 1)/(2*r^4);
% Fs = -3*mu*J2*R_earth^2*(sin(y(2))^2*sin(u)*cos(u))/(r^4);
% Fw = -3*mu*J2*R_earth^2*(sin(u)^2*sin(y(2))*cos(y(2)))/(sin(u)*r^4);
Fr = mu/r^2*(3*J2*(R_earth/r)^2)*P2;
Fs = -mu/r^2*sin(y(2))*sin(u)*J2*(R_earth/r)^2*dP2;
Fw = -mu/r^2*cos(y(2))*J2*(R_earth/r)^2*dP2;

dadt = 2/(n*sqrt(1-y(3)^2))*(y(3)*sin(v)*Fr + p/r*Fs);
dedt = sqrt(1-y(3)^2)/(n*y(1))*(sin(v)*Fr + (cos(v) + (y(3) + cos(v))/(1 + y(3)*cos(v)))*Fs);
didt = r*cos(u)/(n*y(1)^2*sqrt(1-y(3)^2))*Fw;
domegadt = r*sin(u)/(n*y(1)^2*sqrt(1-y(3)^2)*sin(y(2)))*Fw;
dwdt = sqrt(1-y(3)^2)/(n*y(1)*y(3))*(-cos(v)*Fr + sin(v)*(1+r/p)*Fs) - r*cot(y(2))*sin(u)/(sqrt(mu*p))*Fw;
% dvdt = sqrt(mu*p)/r^2 + 1/(y(3)*sqrt(mu*p))*(p*cos(y(6)))*Fr - ((p + r)*sin(y(6)))*Fs;
dndt = -3/2*n/y(1)*dadt;
dMdt = 1/(n*y(1)^2*y(3))*((p*cos(v)-2*y(3)*r)*Fr - (p + r)*sin(v)*Fs) - dndt*(t);

dydt = [dadt; didt; dedt; dwdt; domegadt; dMdt];

end