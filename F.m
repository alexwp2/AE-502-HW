function F = F(z,dt,mu,A,r1mag,r2mag)
    [C,S] = CandS(z);
    F = (y(z,r1mag,r2mag,A)./C).^(1.5).*S + A.*sqrt(y(z,r1mag,r2mag,A)) - sqrt(mu).*dt;
end
