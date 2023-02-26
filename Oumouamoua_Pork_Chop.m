%% AE 502 HW 1 p3 1I/ Oumouamoua Pork Chop Plots
% NOTE: This file takes about 60 seconds to run.
% day 0 is Jan 1, 2017
mu = 1.327E11; % km^3/s^2 [sun]
au_km = 149597870.7; % 1 au to km conversion
day_s = 24*60*60; % 1 day to s conversion

dept_day = 0:1:365;
arr_day = 213:1:760;

% 1I/ Oumouamoua
r1I = [3.515868886595499E-2*au_km, -3.162046390773074*au_km, 4.493983111703389*au_km]; % km
v1I = [-2.317577766980901E-3*au_km/day_s, 9.843360903693031E-3*au_km/day_s, -1.541856855538041E-2*au_km/day_s]; % km/s

% Earth
re = [-1.796136509111975E-1*au_km, 9.667949206859814E-1*au_km, -3.668681017942158E-5*au_km]; % km
ve = [-1.720038360888334E-2*au_km/day_s, -3.211186197806460E-3*au_km/day_s, 7.927736735960840E-7*au_km/day_s]; % km/s

% positions and velocities at departure and arrival points
for i = 1:length(dept_day)
    dt = dept_day(i)*day_s;
    [re_init(i,:),ve_init(i,:)] = Univ_2B_orbit_prop(mu,dt,re,ve);
end
for j = 1:length(arr_day)
    dt = arr_day(j)*day_s;
    [r1I_fin(j,:),v1I_fin(j,:)] = Univ_2B_orbit_prop(mu,dt,r1I,v1I);
end

% direction = 1; % prograde = 1, retrograde = -1
% lambert solutions
for i = 1:length(dept_day)
    for j = 1:length(arr_day)
        if arr_day(j) > dept_day(i)
        dt = (arr_day(j) - dept_day(i))*day_s;
        % prograde
        direction = 1;
        [v1_pro,v2_pro] = Lambert_Solver(mu,dt,re_init(i,:),r1I_fin(j,:),direction);
        % retrograde
        direction = -1;
        [v1_retro,v2_retro] = Lambert_Solver(mu,dt,re_init(i,:),r1I_fin(j,:),direction);
        dv1_pro = v1_pro - ve_init(i,:);
        dv1_retro = v1_retro - ve_init(i,:);
        if dv1_pro <= dv1_retro
            dv1_flyby(i,j) = sqrt(dot(dv1_pro,dv1_pro));
        else
            dv1_flyby(i,j) = sqrt(dot(dv1_retro,dv1_retro));
        end
        dv2_pro = v2_pro - v1I_fin(j,:);
        dv2_retro = v2_retro - v1I_fin(j,:);
        % check which is smaller and use smallest
        if (dv1_pro + dv2_pro) <= (dv1_retro + dv2_retro)
            dv1(i,j) = sqrt(dot(dv1_pro,dv1_pro));
            dv2(i,j) = sqrt(dot(dv2_pro,dv2_pro));
        else
            dv1(i,j) = sqrt(dot(dv1_retro,dv1_retro));
            dv2(i,j) = sqrt(dot(dv2_retro,dv2_retro));
        end
        else
            dv1_flyby(i,j) = 1E5;
            dv1(i,j) = 1E5;
            dv2(i,j) = 1E5;
        end
    end
end

for i = 1:length(dept_day)
    for j = 1:length(arr_day)
        if dv1(i,j)+dv2(i,j) > 50
            dv_tot(i,j) = 50;
        else
            dv_tot(i,j) = dv1(i,j)+dv2(i,j);
        end
    end
end

figure (1)
hold on;
contourf(dept_day,arr_day,transpose(dv_tot),0:.1:50,'LineColor','none');
contourcbar;
title('Delta-V Pork Chop Plot - Rendezvous (dV < 50 km/s)');
xlabel('Departure day (days after Jan 1, 2017)');
ylabel('Arrival day (days after Jan 1, 2017)');

for i = 1:length(dept_day)
    for j = 1:length(arr_day)
        if dv1_flyby(i,j) > 20
            dv_flyby(i,j) = 20;
        else
            dv_flyby(i,j) = dv1_flyby(i,j);
        end
    end
end

figure (2)
hold on;
contourf(dept_day,arr_day,transpose(dv_flyby),0:.1:20,'LineColor','none');
contourcbar;
title('Delta-V Pork Chop Plot - Flyby (dV < 20 km/s)');
xlabel('Departure day (days after Jan 1, 2017)');
ylabel('Arrival day (days after Jan 1, 2017)');

