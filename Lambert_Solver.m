function [v1,v2,h,inc,omega,e,w,theta,a,T] = Lambert_Solver(mu,dt,r1,r2,direction)
%% AE 502 HW 1 p2 - Lambert Solver
% Constants
% mu = 1.327E11; % km^3/s^2 [sun]

% Test values
% mu = 3.986E5; % km^3/s^2 [earth]
% dt = 60*60; % s

% r1 = [5000, 10000, 2100]; % km
% r2 = [-14600, 2500, 7000]; % km

%direction = 1; % prograde = 1, retrograde = -1

r1mag = sqrt(dot(r1,r1));
r2mag = sqrt(dot(r2,r2));

cross12 = cross(r1,r2);
theta = acos(dot(r1,r2)/(r1mag*r2mag))*180/pi;
if direction >= 0
    if cross12(3) <= 0
        theta = 360 - theta;
    end
end
if direction <= 0
    if cross12(3) >= 0
        theta = 360 - theta;
    end
end

A = sin(theta*pi/180)*sqrt(r1mag*r2mag/(1-cos(theta*pi/180)));

z = -20;
while F(z,dt,mu,A,r1mag,r2mag) < 0
    z = z + 0.1;
end

ratio = 1;
while (abs(ratio) > 1E-8)
    ratio = F(z,dt,mu,A,r1mag,r2mag)./dFdz(z,A,r1mag,r2mag);
    z = z - ratio;
end

ymag = y(z,r1mag,r2mag,A);

[C,S] = CandS(z);

f = 1-ymag/r1mag;
g = A*sqrt(ymag/mu);
fdot = sqrt(mu)/(r1mag*r2mag)*sqrt(ymag/C)*(z*S-1);
gdot = 1 - ymag/r2mag;

v1 = 1/g*(r2-f*r1);
v2 = 1/g*(gdot*r2-r1);

[h,inc,omega,e,w,theta,a,T] = orbital_elements(r1,v1,mu);






