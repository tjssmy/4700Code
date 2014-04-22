function GetForces(PhiCutoff,alpha,rm)
global x y nAtoms Fx Fy

for i = 1:nAtoms
    Fx(i) = 0;
    Fy(i) = 0;
    
    for j = 1:nAtoms
        if i == j, continue; end
        dx = x(i)-x(j);
        dy = x(i)-x(j);
        r = sqrt(dx^2 + dy^2);
        
        if r > PhiCutoff, continue, end
        
        [Phi dPhidr] = LJPot(r,alpha,rm);
        
        ang = atan2(dy,dx);
%         dx =  r*cos(ang)
%         dy =  r*sin(ang)
        
        Fx(i) = Fx(i) - dPhidr*cos(ang);
        Fy(i) = Fy(i) - dPhidr*sin(ang);
        

    end
end



end

