function [r,v] = Refined_Orbit_Determination(r2, v2, R_A, Declination, Elevation, Azimuth, MJD)
%% Refined Orbit Determination
% Constants
mu = 3.986E5; % km^3/s^2
day_s = 60*60*24;
Latitude = 40.1164; % deg N
Longitude = -88.2434; % deg E
Height = .233; % km
cluster = 1;

% Initial Orbit Determination
[R1] = Position(Latitude, Longitude, Height, MJD(1 + 3*(cluster - 1))); % km
[R2] = Position(Latitude, Longitude, Height, MJD(2 + 3*(cluster - 1))); % km
[R3] = Position(Latitude, Longitude, Height, MJD(3 + 3*(cluster - 1))); % km
% Direction cosine vector
[rho1] = Direction(R_A(1 + 3*(cluster - 1)), Declination(1 + 3*(cluster - 1)));
[rho2] = Direction(R_A(2 + 3*(cluster - 1)), Declination(2 + 3*(cluster - 1)));
[rho3] = Direction(R_A(3 + 3*(cluster - 1)), Declination(3 + 3*(cluster - 1)));

% Time Intervals
Tau1 = (MJD(1 + 3*(cluster - 1)) - MJD(2 + 3*(cluster - 1)))*day_s;
Tau3 = (MJD(3 + 3*(cluster - 1)) - MJD(2 + 3*(cluster - 1)))*day_s;
Tau = Tau3 - Tau1;

P1 = cross(rho2,rho3);
P2 = cross(rho1,rho3);
P3 = cross(rho1,rho2);

D0 = dot(rho1,P1);
D = [dot(R1,P1), dot(R1,P2), dot(R1,P3); dot(R2,P1), dot(R2,P2), ...
    dot(R2,P3); dot(R3,P1), dot(R3,P2), dot(R3,P3)];

A = 1./D0.*(-D(1,2).*Tau3./Tau + D(2,2) + D(3,2).*Tau1./Tau);
B = 1./(6.*D0).*(D(1,2).*(Tau3^2 - Tau^2).*Tau3./Tau + D(3,2).*(Tau^2 - Tau1^2).*Tau1./Tau);
E = dot(R2,rho2);
R2_square = dot(R2,R2);

a = -(A.^2 + 2.*A.*E + R2_square);
b = -2.*mu.*B.*(A + E);
c = -mu.^2.*B.^2;

Roots = roots([1 0 a 0 0 b 0 0 c]);
x = Roots(Roots >= 0 & imag(Roots) == 0);

rho1_mag = 1/D0*((6*(D(3,1)*Tau1/Tau3 + D(2,1)*Tau/Tau3)*x^3 + mu*D(3,1)*(Tau^2 - Tau1^2)*Tau1/Tau3)/ ...
    (6*x^3 + mu*(Tau^2 - Tau3^2)) - D(1,1));
rho2_mag = A + mu*B/(x^3);
rho3_mag = 1/D0*((6*(D(1,3)*Tau3/Tau1 - D(2,3)*Tau/Tau1)*x^3 + mu*D(1,3)*(Tau^2 - Tau3^2)*Tau3/Tau1)/ ...
    (6*x^3 + mu*(Tau^2 - Tau3^2)) - D(3,3));

% Refined Orbit Determination
rho1_old = rho1_mag; 
rho2_old = rho2_mag; 
rho3_old = rho3_mag;
diff1 = 1; 
diff2 = 1; 
diff3 = 1;

n = 0;
nmax = 1000;
tol = 1.e-8;
% Iterative improvement loop
while ((diff1 > tol) && (diff2 > tol) && (diff3 > tol)) && (n < nmax)
    n = n+1;
    r2_mag = sqrt(dot(r2,r2));
    v2_mag = sqrt(dot(v2,v2));
    alpha = 2/r2_mag - v2_mag^2/mu;
    v2_r = dot(v2,r2)/r2_mag;

    X1 = univ_anomaly(mu,Tau1,r2_mag,v2_r,alpha);
    X3 = univ_anomaly(mu,Tau3,r2_mag,v2_r,alpha);

    [C1,S1] = CandS(alpha*X1^2);
    f1 = 1 - X1^2/r2_mag*C1;
    g1 = Tau1 - X1^3/mu*S1;
    [C3,S3] = CandS(alpha*X3^2);
    f3 = 1 - X3^2/r2_mag*C3;
    g3 = Tau3 - X3^3/mu*S3;

    c1 = g3/(f1*g3 - f3*g1);
    c3 = -g1/(f1*g3 - f3*g1);

    rho1_mag = 1/D0*(-D(1,1) + 1/c1*D(2,1) - c3/c1*D(3,1));
    rho2_mag = 1/D0*(-c1*D(1,2) + D(2,2) - c3*D(3,2));
    rho3_mag = 1/D0*(-c1/c3*D(1,3) + 1/c3*D(2,3) - D(3,3));

    r1 = R1 + rho1_mag*rho1;
    r2 = R2 + rho2_mag*rho2;
    r3 = R3 + rho3_mag*rho3;

    v2 = 1./(f1*g3 - f3*g1).*(-f3.*r1 + f1.*r3);

    diff1 = abs(rho1_mag - rho1_old);
    diff2 = abs(rho2_mag - rho2_old);
    diff3 = abs(rho3_mag - rho3_old);

    rho1_old = rho1_mag; 
    rho2_old = rho2_mag; 
    rho3_old = rho3_mag;

end
r = r2;
v = v2;
end