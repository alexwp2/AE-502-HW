%% Two Body Problem
function [dr_vdt] = twobsolver(t,r_v)
mu = 3.986E5; % km^3/s^2
% Position
r = [r_v(1), r_v(2), r_v(3)];
R = sqrt(dot(r,r));

% Velocity
v = [r_v(4), r_v(5), r_v(6)]';

% Acceleration
a = [-mu*r(1)/R^3,  -mu*r(2)/R^3,  -mu*r(3)/R^3]';

dr_vdt = [v;a];
end