function y = y(z,r1mag,r2mag,A)
    [C,S] = CandS(z);
    y = r1mag + r2mag + A.*(z.*S - 1)./sqrt(C);
end