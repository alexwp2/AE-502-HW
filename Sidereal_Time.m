function [LST] = Sidereal_Time(MJD,Longitude)
JD = MJD + 2400000.5; % Convert MJD to JD
timeDate = datetime(JD,'convertfrom','juliandate','Format','MM/dd/yyyy HH:mm:ss');
[h, m, s] = hms(timeDate);
UT = h + m/60 + s/3600; % UT 

% Local Sidereal Time
T0 = (JD - 2451545)/36525;
GST0 = 100.4606184 + 36000.77004*T0 + 0.000387933*T0^2 - 2.583E-8*T0^3; % Greenwich sidereal time [deg]
GST0 = rem(GST0,360);
GST = GST0 + 360.98564724*UT/24;
LST = GST + Longitude;
end