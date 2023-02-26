%% AE 502 HW 1 p5 I2/ Borisov Plotted Orbit
mu = 1.327E11; % km^3/s^2 [sun]
au_km = 149597870.7; % 1 au to km conversion
day_s = 24*60*60; % 1 day to s conversion
time = 1:1:1855;
dt = time.*day_s;

% 2I/ Borisov
r2I = [7.249472033259724*au_km, 14.61063037906177*au_km, 14.24274452216359*au_km]; % km
v2I = [-8.241709369476881E-3*au_km/day_s, -1.156219024581502E-2*au_km/day_s, -1.317135977481448E-2*au_km/day_s]; % km/s

r2I_fin(1,:) = r2I;
for i = 1:length(time)
    [r2I_fin(i+1,:),~] = Univ_2B_orbit_prop(mu,dt(i),r2I,v2I);
end

dir = 1;
[v1,v2,h,inc,omega,e,w,theta,a,T] = Lambert_Solver(mu,dt(800),r2I,r2I_fin(800,:),dir)

%% Earth Orbit
inc = 0 * pi/180;
OMEGA = 0 * pi/180;
omg = 0 * pi/180;

a = 1.496e+8;
e = 0;
p = a * (1 - e^2);

n = 1;
for t_star = 0:.01:2*pi
    n = n + 1;
    r(n) = p / (1 + e * cos(t_star));
    theta = omg + t_star;
    r_hat(n,1) = t_star;
    r_hat(n,2) = cos(OMEGA) * cos(theta) - sin(OMEGA) * cos(inc) * sin(theta);
    r_hat(n,3) = -1 * cos(OMEGA) * sin(theta) - sin(OMEGA) * cos(inc) * cos(theta);
    r_hat(n,4) = sin(OMEGA) * sin(inc);
    r_x(n) = r(n) * r_hat(n,2);
    r_y(n) = r(n) * r_hat(n,3);
    r_z(n) = r(n) * r_hat(n,4);
end
r_x(1) = r_x(length(r_x));
r_x(2) = r_x(length(r_x));
r_y(1) = r_y(length(r_y));
r_y(2) = r_y(length(r_y));
r_z(1) = r_z(length(r_z));
r_z(2) = r_z(length(r_z));

figure (1)
hold on;
axis equal;
plot3(0,0,0,'r.','LineWidth',30);
plot3(r2I_fin(:,1),r2I_fin(:,2),r2I_fin(:,3),'b-');
plot3(r_x,r_y,r_z,'g-');

