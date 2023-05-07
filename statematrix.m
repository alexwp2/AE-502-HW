function [dphidt] = statematrix(t, r_v, F)
phi = zeros(6);
for i = 1:6
    phi(:,i) = r_v(1+6*(i-1):6+6*(i-1));
end
dphi = F*phi;

dphidt = [dphi(:,1); dphi(:,2); dphi(:,3); dphi(:,4); dphi(:,5); dphi(:,6)];
end