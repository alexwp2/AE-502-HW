function [X] = univ_anomaly(mu,dt,r_0_mag,vr0,alpha)
    i = 1;
    Xi(i) = sqrt(mu)*abs(alpha)*dt;
    ratio = 1;
    while abs(ratio) > 1E-8
        z = alpha*Xi(i)^2;
        [C,S] = CandS(z);     
        fX(i) = r_0_mag*vr0/sqrt(mu)*Xi(i)^2*C + (1-alpha*r_0_mag)*Xi(i)^3*S + r_0_mag*Xi(i) - sqrt(mu)*dt;
        dfX(i) = r_0_mag*vr0/sqrt(mu)*Xi(i)*(1-alpha*Xi(i)^2*S) + (1-alpha*r_0_mag)*Xi(i)^2*C + r_0_mag;
        ratio = fX(i)/dfX(i);
        if abs(ratio) > 1E-8
            Xi(i+1) = Xi(i) - ratio;
        end
        i = i + 1;
    end
    X = Xi(length(Xi));

end