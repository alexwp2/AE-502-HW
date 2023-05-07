function [r_EKF,v_EKF] = ExtKalman(R_A, Declination, Elevation, Azimuth, MJD)
%% Extended Kalman Filter
% Constants
mu = 3.986E5; % km^3/s^2
day_s = 60*60*24;
Latitude = 40.1164; % deg N
Longitude = -88.2434; % deg E
Height = .233; % km
error_range = .001; % km
error_Azimuth = 0.000277778*pi/180; % 1 arc sec in rad
error_Elevation = 0.000277778*pi/180; % 1 arc sec in rad

% Average r and v of the three clusters
% cluster = 1;
% [r1, v1] = Initial_Orbit_Determination(R_A, Declination, Elevation, Azimuth, MJD, cluster);
% cluster = 2;
% [r2, v2] = Initial_Orbit_Determination(R_A, Declination, Elevation, Azimuth, MJD, cluster);
% cluster = 3;
% [r3, v3] = Initial_Orbit_Determination(R_A, Declination, Elevation, Azimuth, MJD, cluster);
% 
% r = [r1(1) + r2(1) + r3(1), r1(2) + r2(2) + r3(2), r1(3) + r2(3) + r3(3)]./3;
% v = [v1(1) + v2(1) + v3(1), v1(2) + v2(2) + v3(2), v1(3) + v2(3) + v3(3)]./3;
% r_v = [r(1), r(2), r(3), v(1), v(2), v(3)];

[r2, v2] = Initial_Orbit_Determination(R_A, Declination, Elevation, Azimuth, MJD, 1);
[r,v] = Refined_Orbit_Determination(r2, v2, R_A, Declination, Elevation, Azimuth, MJD);
x_state = [r(1), r(2), r(3), v(1), v(2), v(3)];

% Initial epoch
t0 = (MJD(2) - MJD(1))*day_s;
tspan0 = [0, t0];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~, x_ref] = ode45(@twobsolver, tspan0, x_state, options);

x_init = x_ref(end,:);
Weights = [(1/error_range)^2, 0, 0; 0, (1/error_Azimuth)^2, 0; 0, 0, (1/error_Elevation)^2];
N = length(MJD)/9; % Number of observations;

% Initial Covariance Matrix
P = zeros(6);
for i=1:3
    P(i,i)=1e8;
end
for i=4:6
    P(i,i)=1e3;
end

x_store = transpose(x_init);
for i = 1:N
    x_old = x_store;
    phi = eye(6,6);
    if i == 1
        t = -t0;
    else
        t = (MJD(i)-MJD(i-1))*24*60*60;
    end
    % Propagate
    tspan = [0 t];
    [~, x_med] = ode45(@twobsolver, tspan, x_old, options);
    x_store = transpose(x_med(end,:));

    r_mag = sqrt(dot(x_store(1:3),x_store(1:3)));
    F = zeros(6);
    F(1:3,4:6) = eye(3);
    for i = 1:3
        for j = 1:3
            if i == j
                F(i+3,j) = (-mu/r_mag^3)+((3*mu*x_store(i)^2)/r_mag^5);
            end
            F(i+3,j) = (3*mu*x_store(i)*x_store(j))/r_mag^5;
        end
    end
    [~,phi_store] = ode45(@statematrix, tspan, phi, options, F);
    phi_new = phi_store(end,:);
    for j = 1:6
        phi(:,j) = phi_new(1 + 6*(j - 1):6 + 6*(j - 1));
    end
    P = phi*P*transpose(phi);
    
    m = 0;
    while i - 3 > 0
        m = m + 1;
    end
    
    [R1] = Position(Latitude, Longitude, Height, MJD(1 + 3*(m))); % km
    [R2] = Position(Latitude, Longitude, Height, MJD(2 + 3*(m))); % km
    [R3] = Position(Latitude, Longitude, Height, MJD(3 + 3*(m))); % km

    % Direction cosine vector
    [rho1] = Direction(R_A(1 + 3*(m)), Declination(1 + 3*(m)));
    [rho2] = Direction(R_A(2 + 3*(m)), Declination(2 + 3*(m)));
    [rho3] = Direction(R_A(3 + 3*(m)), Declination(3 + 3*(m)));

    Tau1 = (MJD(1 + 3*(m)) - MJD(2 + 3*(m)))*day_s;
    Tau3 = (MJD(3 + 3*(m)) - MJD(2 + 3*(m)))*day_s;
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

    r_site = rho2_mag*rho2;

    range_predict = norm(r_site);
    rho_predict = x_store(1:3) - r_site(1:3);
%     range_predict = sqrt(dot(rho_predict,rho_predict));
%     unit_rho = rho_predict/range_predict;
    unit_rho = -rho2;
    %% Debugged to here
    Declination_predict = asind(unit_rho(2));
    R_A_predict = acosd(unit_rho(1)/cosd(Declination_predict));
    if rho2(2) <= 0
        R_A_predict = 360 - R_A_predict;
    end
    [LST] = Sidereal_Time(MJD(i),Longitude);
    Sid_Time = mod(LST - R_A_predict, 360);
    Elevation_predict = asind(sind(Latitude).*sind(Declination_predict)+cosd(Latitude).*cosd(Declination_predict).*cosd(Sid_Time));
    Azimuth_predict = mod(atan2(-sind(Sid_Time).*cosd(Declination_predict)./cosd(Elevation_predict),...
        (sind(Declination_predict)-sind(Elevation_predict).*sind(Latitude))./(cosd(Elevation_predict).*cosd(Latitude))).*(180/pi),360);

    % Observation Matrix H
    Rot3 = [cos(-deg2rad(LST)), sin(-deg2rad(LST)), 0; -sin(-deg2rad(LST)), cos(-deg2rad(LST)), 0; 0, 0, 1];
    Rot2 = [cos(-((pi/2)-deg2rad(Longitude))), 0, -sin(-((pi/2)-deg2rad(Longitude))); 0, 1, 0; sin(-((pi/2)-deg2rad(Longitude))), 0, cos(-((pi/2)-deg2rad(Longitude)))];
    Rot_matrix = Rot3*Rot2;

    rho_rot_predict = transpose(Rot_matrix).*rho_predict;
    rho_rot_predict_mag = norm(rho_rot_predict);
    rho_rot_12 = norm(rho_rot_predict(1:2));

    obs_rho = [rho_rot_predict(1)/rho_rot_predict_mag, rho_rot_predict(2)/rho_rot_predict_mag, rho_rot_predict(3)/rho_rot_predict_mag;
        -rho_rot_predict(2)/rho_rot_12^2, rho_rot_predict(1)/rho_rot_12^2, 0; rho_rot_predict(1)*rho_rot_predict(3)/...
        ((rho_rot_predict_mag^2)*rho_rot_12), rho_rot_predict(2)*rho_rot_predict(3)/((rho_rot_predict_mag^2)*rho_rot_12), -rho_rot_12./(rho_rot_predict_mag^2)];
    
    H = [obs_rho*transpose(Rot_matrix), zeros(3)];

    % residual
    observed = horzcat(horzcat(rho2_mag,Azimuth),Elevation); 
    rad_obs = deg2rad(observed(:,2:3));
    residual = [observed(1,1) rad_obs(1,:)]' - transpose([range_predict deg2rad(Azimuth_predict) deg2rad(Elevation_predict)]);
    K = P*transpose(H)*pinv(H*P*transpose(H) + inv(Weights));

    delx = K*(residual);
    x_store = x_store + delx;
    P = (eye(6) - K*H)*P;
end

tspanzero = [0, (MJD(end) - MJD(1))*day_s];
[~, xzero] = ode45(@twobsolver, tspanzero, x_store, options);
EKF = transpose(xzero(end,:));

r_EKF = [EKF(1), EKF(2), EKF(3)];
v_EKF = [EKF(4), EKF(5), EKF(6)];
% r_EKF = x_state(1:3);
% v_EKF = x_state(4:6);


end