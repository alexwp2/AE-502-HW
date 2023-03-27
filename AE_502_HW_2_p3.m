%% AE 502 HW 2 problem 3
% Constants
mu = 3.986E5; % km^3/s^2 [earth]
R_earth = 6370; % km
J2 = .00108; % J2 perturbation
day = 60*60*24; % s

% Test Initial Conditions
% tspan = 100*day/24;
% a_0 = 8309; % km
% i_0 = 0.4887; % rad
% e_0 = 0.1963;
% w_0 = 30*pi/180; % rad
% omega_0 = 45*pi/180; % rad
% M_0 = 40*pi/180; % rad



% Initial Conditions
a_0 = 26600; % km
i_0 = 1.10654; % rad
e_0 = 0.74;
w_0 = 5*pi/180; % drad
omega_0 = 90*pi/180; % rad
M_0 = 10*pi/180; % rad
Period = 2*pi*sqrt(a_0^3/mu); % sec
tspan = 20*day;


% E_0 = M_0 + e_0*sin(M_0) + e_0^2/2*sin(2*M_0) + e_0^3/8*(3*sin(3*M_0) - sin(M_0)) + e_0^4/6*(2*sin(4*M_0) - sin(2*M_0));
% v_0 = 2*atan(sqrt((1+e_0)/(1-e_0))*tan(E_0/2));
r = a_0.*(1 - e_0^2)./(1 + e_0.*cos(0:2*pi/day:101*2*pi));

y_0 = [a_0, i_0, e_0, w_0, omega_0, M_0];
% Solution
opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t,y] = ode45(@(t,y) GaussPlanetaryEqs(t,y,y_0,mu,J2,R_earth,r,Period), [0:60:tspan], y_0, opts);
len = length(y(:,1));
% opts = odeset('RelTol',1e-4,'AbsTol',1e-4);
% [t,y] = ode45(@(t,y) LagrangePlanetaryEqs(t,y,y_0,mu,J2,R_earth,r,Period), [0:60:tspan], y_0, opts);

figure (4)
hold on;
plot(t./day,y(:,4)*180/pi);
title('w Perturbation');
xlabel('Time (days)');
ylabel('w (degrees)');

figure (5)
hold on;
plot(t./day,y(:,5)*180/pi);
title('Omega Perturbation');
xlabel('Time (days)');
ylabel('Omega (degrees)');

for x = 1:20
    for j = 1:(len-1)/20
    y(x*(len-1)/20 + j,1) = y(j,1);
    y(x*(len-1)/20 + j,2) = y(j,2);
    y(x*(len-1)/20 + j,3) = y(j,3);
    t(x*(len-1)/20 + j) = t(j) + day*x;
    end
end


figure (1)
hold on;
plot(t./day,y(:,1));
xlim([0,20]);
title('Semi-major Axis Perturbation');
xlabel('Time (days)');
ylabel('a (km)');

figure (2)
hold on;
plot(t./day,y(:,2)*180/pi);
xlim([0,20]);
title('Inclination Perturbation');
xlabel('Time (days)');
ylabel('Inclination (degrees)');

figure (3)
hold on;
plot(t./day,y(:,3));
xlim([0,20]);
title('Eccentricity Perturbation');
xlabel('Time (days)');
ylabel('e');





