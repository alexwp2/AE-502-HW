function dF = dFdz(z,A,r1mag,r2mag)
    [C,S] = CandS(z);
    if z == 0
        dF = sqrt(2)./40.*y(z,r1mag,r2mag,A).^(1.5) + A./8.*(sqrty(z,r1mag,r2mag,A)) ...
        + A.*sqrt(1./(2.*y(z,r1mag,r2mag,A)));
    else
        dF = (y(z,r1mag,r2mag,A)./C).^(1.5).*(1./(2.*z).*(C - 3.*S./(2.*C)) + 3.*S.^2./(4.*C)) ...
        + A./8.*(3.*S./C.*sqrt(y(z,r1mag,r2mag,A)) + A.*sqrt(C./y(z,r1mag,r2mag,A)));
    end
end