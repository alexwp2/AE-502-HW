%% AE 502 HW 4 Main
clear;
% Constants
mu = 3.986E5; % km^3/s^2
day_s = 60*60*24;

% Initial Conditions
R_A = [265.3869518, 286.0157158, 306.8358871, 113.1072289, 44.4521793, 27.5410256, 128.5555539, 116.9894450, 102.4669647]; % [deg] Right Ascension
Declination = [38.4539569, 47.7809854, 50.4514381, 76.1961378, 71.6883545, 61.8733178, 21.4050426, 32.8319850, 41.0669857]; % [deg] Declination
Elevation = [19.3547079, 14.2703198, 8.1356679, 40.4709053, 23.6357032, 12.4517325, 10.1299902, 9.4168870, 6.8773277]; % [deg] Elevation
Azimuth = [55.5481428, 38.2354388, 25.9457472, 341.8114490, 350.7827597, 355.0101055, 289.4238109, 305.2967884, 319.7464073]; % [deg] Azimuth
MJD = [60055.1340277, 60055.1354166, 60055.1368055, 60055.2104166, 60055.2118055, 60055.2131944, 60055.2854166, 60055.2868055, 60055.2881944]; % Modified Julian Date

% Initial Conditions (Cluster 1)
% R_A = [265.3869518, 286.0157158, 306.8358871]; % [deg] Right Ascension
% Declination = [38.4539569, 47.7809854, 50.4514381]; % [deg] Declination
% Elevation = [19.3547079, 14.2703198, 8.1356679]; % [deg] Elevation
% Azimuth = [55.5481428, 38.2354388, 25.9457472]; % [deg] Azimuth
% MJD = [60055.1340277, 60055.1354166, 60055.1368055]; % Modified Julian Date

% Initial Conditions (Cluster 2)
% R_A = [113.1072289, 44.4521793, 27.5410256]; % [deg] Right Ascension
% Declination = [76.1961378, 71.6883545, 61.8733178]; % [deg] Declination
% Elevation = [40.4709053, 23.6357032, 12.4517325]; % [deg] Elevation
% Azimuth = [341.8114490, 350.7827597, 355.0101055]; % [deg] Azimuth
% MJD = [60055.2104166, 60055.2118055, 60055.2131944]; % Modified Julian Date

% % Initial Conditions (Cluster 3)
% R_A = [128.5555539, 116.9894450, 102.4669647]; % [deg] Right Ascension
% Declination = [21.4050426, 32.8319850, 41.0669857]; % [deg] Declination
% Elevation = [10.1299902, 9.4168870, 6.8773277]; % [deg] Elevation
% Azimuth = [289.4238109, 305.2967884, 319.7464073]; % [deg] Azimuth
% MJD = [60055.2854166, 60055.2868055, 60055.2881944]; % Modified Julian Date


% Initial Orbit r2, v2, and orbital elements
cluster = 1;
[r2, v2] = Initial_Orbit_Determination(R_A, Declination, Elevation, Azimuth, MJD, cluster);
r2_mag = sqrt(dot(r2,r2));

[~,inc,omega,emag,w,theta,a,T] = orbital_elements(r2,v2,mu);

% Refined Orbit r, v, and orbital elements
[r,v] = Refined_Orbit_Determination(r2, v2, R_A, Declination, Elevation, Azimuth, MJD);
r_mag = sqrt(dot(r,r));
v_mag = sqrt(dot(v,v));
[~,inc_ref,omega_ref,emag_ref,w_ref,theta_ref,a_ref,T_ref] = orbital_elements(r,v,mu);

% Extended Kalman Filter
[r_EKF,v_EKF] = ExtKalman(R_A, Declination, Elevation, Azimuth, MJD);

% Determine Satellite
% Propagate orbit out to 2023-05-05 21:00:00.000 and compare with available
% data on Privateer Wayfinder and available TLE
% 
% Determinte ONEWEB-0XXX orbit rs1, vs1, and orbital elements

% ONEWEB-0113 @ 12:00:00
% rs1 = [-6576.109475227599 -866.8765902693341 3681.229281133656];
% rs2 = [-6109.842677568626 -842.1484444906677 4414.867374882107];
% ONEWEB-0350 2023-05-05 21:00:00.000
% rs1 = [-7498.095919428197 -850.5477259191871 812.6500277355466];
% rs2 = [-7351.946257909613 -866.193002585712 1669.5581394041212];
% ONEWEB-0178 2023-05-05 21:00:00.000
% rs1 = [-6752.823865948233 -860.4758581845351 3348.2543326024547];
% rs2 = [-6323.158456597478 -841.511743849525 4104.0878842951715];
% ONEWEB-0336 2023-05-06 09:43:12.000
% rs1 = [-3713.018062044482 -641.1585622367758 6578.0041934811825];
% rs2 = [-2936.969130866953 -572.2012345530869 6964.398085879295];
% ONEWEB-0112 2023-05-06 13:47:12.000
rs1 = [4473.104430231738 250.33346268666895 6115.981692110349];
rs2 = [5139.25769584076 341.6793777187763 5564.378118448044];

% ONEWEB-0169 2023-05-06 19:06:15.000
% rs1 = [1729.4368967497242 -89.86914386910193 7378.744996964597];
% rs2 = [2558.5564378675463 6.801946507147048 7134.81492785654];
% ONEWEB-0345 2023-05-06 13:50:13.000
% rs1 = [4452.947601678044 247.594567512955 6131.407326037197];
% rs2 = [5120.9421440881415 339.0345620751399 5582.1169748695465];

Endtime = datetime('2023-05-06 13:47:12.000', 'InputFormat', 'uuuu-MM-dd HH:mm:ss.sss');
JD2 = juliandate(Endtime);
MJD2 = JD2 - 2400000.5;
Runtime = (MJD2 - MJD(2))*day_s;
diff = 120; % s

[vs1,~,~,inc_s,omega_s,e_s,w_s,theta_s,a_s,T_s] = Lambert_Solver(mu,diff,rs1,rs2,1);
times = 1:1:T_s;
dts = times;
rs_prop(1,:) = rs1;
for m = 1:length(times)
    [rs_prop(m+1,:),~] = Univ_2B_orbit_prop(mu,dts(m),rs1,vs1);
end

% Plot inital orbit determination
time = 1:1:T;
dt = time;
r_prop(1,:) = r2;
for m = 1:length(time)
    [r_prop(m+1,:),~] = Univ_2B_orbit_prop(mu,dt(m),r2,v2);
end

% Plot refined orbit determination
time2 = 1:1:T_ref;
dt2 = time2;
r_prop_ref(1,:) = r;
for m = 1:length(time2)
    [r_prop_ref(m+1,:),~] = Univ_2B_orbit_prop(mu,dt2(m),r,v);
end

% Plot Extended Kalman Filter orbit determination
[~,inc_EKF,omega_EKF,emag_EKF,w_EKF,theta_EKF,a_EKF,T_EKF] = orbital_elements(r_EKF,v_EKF,mu)
time_EKF = 1:1:T_EKF;
dt_EKF = time_EKF;
r_prop_EKF(1,:) = r_EKF;
for m = 1:length(time_EKF)
    [r_prop_EKF(m+1,:),~] = Univ_2B_orbit_prop(mu,dt_EKF(m),r_EKF,v_EKF);
end

% Propagate Given Satellite to 2023-05-06 13:47:12.000
r_prop_end(1,:) = r;
for i = 1:10
    time = Runtime/10*(i-1) + 1:1:Runtime/10*i;
    dt = time;
    for m = 1:length(time)
        [r_prop_end(m+1,:),~] = Univ_2B_orbit_prop(mu,dt(m),r,v);
    end
    r_prop_end(1,:) = r_prop_end(length(r_prop_end),:);
end
r_end = r_prop_end(length(r_prop_end),:);

% Reverse propagation of ONEWEB-0XXX to initial given satellite data time
r_rev(1,:) = rs1;
time = 1:-5:-Runtime + 1*day_s/24*5.013 + 120;
dt = time;
for m = 1:length(time)
    [r_rev(1,:),~] = Univ_2B_orbit_prop(mu,dt(m),rs1,vs1);
end
r_final = r_rev(1,:);


%% Plots
[x,y,z] = sphere;
x = x*6371;
y = y*6371;
z = z*6371;

ra1 = [-4.9539E3, -1.0292E3, 5.6450E3];
ra3 = [-3.5504E3, -0.8439E3, 6.6357E3];

figure (1)
hold on;
axis equal;
plot3(r_prop(:,1), r_prop(:,2), r_prop(:,3), 'b');  % Given satellite initial determined orbit propagation
plot3(r_prop_ref(:,1), r_prop_ref(:,2), r_prop_ref(:,3), 'r'); % Given satellite refined determined orbit propagation
% plot3(r(1), r(2), r(3),'g.'); % Given satellite initial position r2
% plot3(ra1(1), ra1(2), ra1(3),'gx'); % Given satellite initial position r1
% plot3(ra3(1), ra3(2), ra3(3),'rx'); % Given satellite initial position r3
% plot3(r_end(1), r_end(2), r_end(3), 'k.'); % Given satellite final position
% plot3(rs1(1), rs1(2), rs1(3), 'c.'); % Test satellite final position
% plot3(real(r_final(1)), real(r_final(2)), real(r_final(3)), 'o'); % Test satellite reverse prop
% plot3(real(rs_prop(:,1)), real(rs_prop(:,2)), real(rs_prop(:,3)), 'm');
% plot3(r_EKF(1), r_EKF(2), r_EKF(3),'g.'); % Given satellite initial position
plot3(r_prop_EKF(:,1), r_prop_EKF(:,2), r_prop_EKF(:,3), 'm'); % Given satellite EKF propagation
surf(x,y,z);


