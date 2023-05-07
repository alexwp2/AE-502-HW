function [rho] = Direction(R_A, Declination)
% topocentric equitorial coordinate system
rho = [cos(Declination*pi/180)*cos(R_A*pi/180), cos(Declination*pi/180)*sin(R_A*pi/180), sin(Declination*pi/180)];

% topocentric horizon coordinate system
% rho_horiz = [cos(Elevation*pi/180)*sin(Azimuth*pi/180), cos(Elevation*pi/180)*cos(Azimuth*pi/180), sin(Elevation*pi/180)];
% rho_i = rho_horiz(1).*[-sin(Long*pi/180), cos(Long*pi/180), 0];
% rho_j = rho_horiz(2).*[-sin(Lat*pi/180)*cos(Long*pi/180), -sin(Lat*pi/180)*sin(Long*pi/180), cos(Lat*pi/180)];
% rho_k = rho_horiz(3).*[cos(Lat*pi/180)*cos(Long*pi/180), cos(Lat*pi/180)*sin(Long*pi/180), sin(Lat*pi/180)];
% rho = rho_i + rho_j + rho_k;
end