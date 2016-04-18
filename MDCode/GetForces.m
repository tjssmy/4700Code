function GetForces(PhiCutoff,Epsilon,sigma)
global x y nAtoms Fx Fy Phi

for i = 1:nAtoms
    dx = x(i) - x;
    dy = y(i) - y;
    r = sqrt(dx.^2 + dy.^2);
    
    withinPhi = (r ~= 0) & (r <= PhiCutoff);
    r = r(withinPhi);
    dx = dx(withinPhi);
    dy = dy(withinPhi);
    
    [aPhi, dPhidr] = LJPot(r, Epsilon, sigma);
    ang = atan2(dy, dx);
    dFx = -dPhidr .* cos(ang);
    dFy = -dPhidr .* sin(ang);
    
    Phi(i) = sum(aPhi);
    Fx(i) = sum(dFx);
    Fy(i) = sum(dFy);
end

hold off

end
