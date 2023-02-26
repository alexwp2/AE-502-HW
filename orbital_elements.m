function [h,inc,omega,emag,w,theta,a,T] = orbital_elements(r,v,mu)
    rmag = sqrt(dot(r,r));
    vmag = sqrt(dot(v,v));
    vr = dot(r,v)/rmag;
    h = cross(r,v); 
    hmag = sqrt(dot(h,h)); % specific angular momentum
    inc = acos(h(3)/hmag)*180/pi; % inclination
    N = cross([0,0,1],h);
    Nmag = sqrt(dot(N,N));
    omega = acos(N(3)/Nmag)*180/pi; % RA of ascending node
    if N(2) < 0
        omega = 360-omega;
    end
    e = 1/mu*((vmag^2-mu/rmag).*r - rmag.*vr.*v);
    emag = sqrt(dot(e,e)); % eccentricity
    w = acos(dot(N,e)/(Nmag*emag))*180/pi; % argument of periapsis
    if e(3) < 0
        w = 360-w;
    end
    theta = acos(dot(e,r)/(emag*rmag))*180/pi; % true anomaly
    if vr < 0
        theta = 360-theta;
    end
    rp = hmag^2/mu/(1+emag);
    ra = hmag^2/mu/(1-emag);
    a = .5*(rp+ra); % semimajor axis
    T = 2*pi/sqrt(mu)*a^(3/2); % period
end