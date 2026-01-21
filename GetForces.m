function GetForces(PhiCutoff,Epsilon,sigma)
global x y nAtoms Fx Fy Phi
global step applyForce LAtoms WAtoms


squish = false;             %%initially squish is set to false
if applyForce               %%if apply force is true, force will initally be applied - tried to see how interactions are after large force is applied inwards
    if step == 0            %%only applies the force initially (time 0)
        squish = true;             %%will apply force fo 1e-9 to top and bottom electrons (facing inwards) -  atomic array will fall apart
    end
end

% n = 1;
for i = 1:nAtoms            %%cycles through for each atoms
    Fx(i) = 0;              %%initially set everything for that atom to zero
    Fy(i) = 0;
    Phi(i) = 0;
%     if i == 10
%         plot(x(i),y(i),'o','markers',48);
%         hold on
%     end

    for j = 1:nAtoms            %cycles through atoms to calcuate force on the number i atoms
        if i == j, continue; end        %%if i = j it will pass to the next iteration of the for loop - can't ahve a force from itself
        dx = x(i) - x(j);
        dy = y(i) - y(j);
        r = sqrt(dx^2 + dy^2);

        if r > PhiCutoff, continue, end         %if atom is too far away, it does not have an effect on the $i atom

        [aPhi dPhidr] = LJPot(r, Epsilon, sigma);

        ang = atan2(dy, dx);
%         dx =  r*cos(ang)
%         dy =  r*sin(ang)

        dFx = - dPhidr * cos(ang);
        dFy = - dPhidr * sin(ang);

%         if i == 10
%             plot(x(j),y(j),'ro','markers',48);
%         end

        Phi(i) = Phi(i) + aPhi;
        Fx(i) = Fx(i) + dFx;
        Fy(i) = Fy(i) + dFy;   
       
        if squish                   %%apply force to top and bottom layers (force inwards) - ignores the calculated forces
            if i > (LAtoms*WAtoms - LAtoms)
                Fy(i) = -1e-9; 
            elseif i < (LAtoms+1)
                Fy(i) = 1e-9;
            end
    end
end

hold off

end

