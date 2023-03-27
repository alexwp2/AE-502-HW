%% AE 502 HW 2 problem 1
% Constants
mu = 3.986E5; % km^3/s^2 [earth]

R_earth = 6370; % km
J2 = .00108; % J2 perturbation
day = 60*60*24; % s

alt_min = 600; % km



% Solution
Period = day/3; % s
n = 360/Period; % rad/s

a = (Period/(2*pi)*sqrt(mu))^(2/3); % km
alt_max = a; % km
r_p = (alt_min + R_earth):alt_max;
i = acos(sqrt(0.2))*180/pi;

for j = 1:length(r_p)
    e(j) = 1 - (r_p(j)/a);
    arg_asc(j) = -3/2*n*J2*(R_earth/a)^(2)*cos(i*pi/180)/(1-e(j)^(2))^(2);
    arg_asc_deg(j) = arg_asc(j)*day;
end

fprintf('Final Orbital Elements:\n');
fprintf('Semi-major axis (km) = ')
a 
fprintf('\nEccentricity = ')
ecc = e(1)
fprintf('\nInclination (degrees) = ')
i
fprintf('Lowest drift rate (deg/day) = ')
arg = arg_asc_deg(1)


