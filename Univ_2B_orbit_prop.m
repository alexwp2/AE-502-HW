function [r,v] = Univ_2B_orbit_prop(mu,dt,r_0,v_0)
%% AE 502 HW 1 p1 - Universal Variable two-body orbit propagator
% Constants
% mu = 1.327E11; % km^3/s^2 [sun]

% Test values
% mu = 3.986E5; % km^3/s^2 [earth]
% dt = 60*60; % s [delta t] 

% r_0 = [7000, -12124, 0]; % km
% v_0 = [2.6679, 4.621, 0]; % km/s

r_0_mag = sqrt(dot(r_0,r_0));
v_0_mag = sqrt(dot(v_0,v_0));

vr0 = dot(r_0,v_0)/r_0_mag;
alpha = 2/r_0_mag - (v_0_mag)^2/mu;

X = univ_anomaly(mu,dt,r_0_mag,vr0,alpha);

z = alpha*X^2;

[C,S] = CandS(z);

f = 1 - X^2/r_0_mag*C;
g = dt - 1/sqrt(mu)*X^3*S;

r = f*r_0 + g*v_0;
r_mag = sqrt(dot(r,r));

fdot = sqrt(mu)/(r_mag*r_0_mag)*(alpha*X^3*S-X);
gdot = 1-X^2/r_mag*C;

v = fdot*r_0 + gdot*v_0;
v_mag = sqrt(dot(v,v));
