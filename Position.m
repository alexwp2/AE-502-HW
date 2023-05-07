function [R] = Position(Latitude, Longitude, H, MJD)
%% Position on Earth at measured time
% J2000 epoch
% Constants
mu = 3.986E5; % km^3/s^2
f = 0.00335; % Earth Oblateness
R_eq = 6378.137; % km
% Latitude = 40.1164; % deg N
% Longitude = -88.2434; % deg E
% H = .233; % km
[LST] = Sidereal_Time(MJD,Longitude);

% Calculate R position
R_psi = R_eq/(sqrt(1 - (2*f - f^2)*sin(Latitude*pi/180)^2));
R_c = R_psi + H;
R_s = (1 - f)^2*R_psi + H;
R = [R_c*cos(Latitude*pi/180)*cos(LST*pi/180), R_c*cos(Latitude*pi/180)*sin(LST*pi/180), R_s*sin(Latitude*pi/180)]; % km
% R_mag = sqrt(dot(R,R));

end