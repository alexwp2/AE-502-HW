function [C,S] = CandS(z)
    if z > 0
        C = (1-cos(sqrt(z)))/z;
    end
    if z < 0
        C = (cosh(sqrt(-z))-1)/(-z);
    end
    if z == 0
        C = 1/2;
    end
    if z > 0
        S = (sqrt(z)-sin(sqrt(z)))/sqrt(z)^3;
    end
    if z < 0
        S = (sinh(sqrt(-z))-sqrt(-z))/sqrt(-z)^3;
    end
    if z == 0
        S = 1/6;
    end  
end