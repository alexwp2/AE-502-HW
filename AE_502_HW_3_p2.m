%% AE 502 HW 3 problem 2
a = 1; % semi-major axis
e = 0.5; % eccentricity
i = 45*pi/180; % inclination [deg] to radians
w = .01;
time = 0:1:100; % time units
Period = 2*pi*a^(3/2); % period
n = 2*pi/Period;

omg = [0, w, 0; -w, 0, 0; 0, 0, 0];

for j = 1:length(time)
    r(:,j) = exp(-time(j)*omg)*[-(a*(1-e) - time(j)/(n*a^2)^3 - (2*n*a^2 +...
        3*n*a^2*sqrt(1-e^2) + 3*n*a^2*sqrt(1-e^2)*cos(i))*time(j)/(n*a^2)^4);...
        (0 - time(j)/(n*a^2)^3); (0 + 1/2 - time(j)/(n*a^2)^3)];
    rmag(j) = sqrt(dot(r(:,j),r(:,j)));
    R(:,j) = r(:,j)./rmag(j);
end

%% Approximate Unperturbed Orbit (for comparison)
inc = 45 * pi/180;
OMEGA = 0 * pi/180;
omg = 0 * pi/180;
p = a * (1 - e^2);

n = 1;
for t_star = 0:.01:2*pi
    n = n + 1;
    r_mag(n) = p / (1 + e * cos(t_star));
    theta = omg + t_star;
    r_hat(n,1) = t_star;
    r_hat(n,2) = cos(OMEGA) * cos(theta) - sin(OMEGA) * cos(inc) * sin(theta);
    r_hat(n,3) = -1 * cos(OMEGA) * sin(theta) - sin(OMEGA) * cos(inc) * cos(theta);
%     r_hat(n,4) = sin(OMEGA) * sin(inc);
    r_x(n) = r_mag(n) * r_hat(n,2);
    r_y(n) = r_mag(n) * r_hat(n,3);
%     r_z(n) = r_mag(n) * r_hat(n,4);
end
r_x(1) = r_x(length(r_x));
r_x(2) = r_x(length(r_x));
r_y(1) = r_y(length(r_y));
r_y(2) = r_y(length(r_y));
% r_z(1) = r_z(length(r_z));
% r_z(2) = r_z(length(r_z));

%% Plots
figure (1)
hold on;
plot3(R(1,:),R(2,:),R(3,:)); % Hamiltionian EOMs
plot3(r_x,r_y,r_z); % Approximate Unperturbed Orbit
plot3(0,0,0,'g.','LineWidth',30); % Origin


