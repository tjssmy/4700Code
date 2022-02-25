function AddHCPAtomicBlob(LAtoms, X0, Y0, VX0, VY0, RotAng, InitDist, Temp, Type)
global C
global x y AtomSpacing
global nAtoms % MinX MaxX MinY MaxY
global AtomType Vx Vy Mass0 Mass1

if Type == 0
    Mass = Mass0;
else
    Mass = Mass1;
end

% L = ((LAtoms-1)+0.5)*AtomSpacing;
% W = (LAtoms-1)*sqrt(3)/2*AtomSpacing;

p = 1;

for j=1:LAtoms + 1
    for i = 1:2*LAtoms - 1 - j
        y(end + 1) = (sqrt(3) * ( j - 1) * AtomSpacing) / 2.0;
        x(end + 1) = (-(2 * LAtoms - j - 2) * AtomSpacing) / 2.0 + i * AtomSpacing;
        p = p + 1;
        if j ~= 1
            y(end + 1) = -y(end);
            x(end + 1) = x(end);
            p = p + 1;
        end
    end
end

x = x - (max(x) + min(x)) / 2;

numAtoms = p-1;

x(nAtoms + 1:nAtoms + numAtoms) = x(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms)-0.5)*AtomSpacing*InitDist + X0*AtomSpacing;
y(nAtoms + 1:nAtoms + numAtoms) = y(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms)-0.5)*AtomSpacing*InitDist + Y0*AtomSpacing;

AtomType(nAtoms + 1:nAtoms + numAtoms) = Type;

if Temp == 0
    Vx(nAtoms + 1:nAtoms + numAtoms) = 0;
    Vy(nAtoms + 1:nAtoms + numAtoms) = 0;
else
    std0 = sqrt(C.kb*Temp/Mass);

    Vx(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
    Vy(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
end

Vx(nAtoms + 1:nAtoms + numAtoms) = Vx(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vx(nAtoms + 1:nAtoms + numAtoms)) + VX0;
Vy(nAtoms + 1:nAtoms + numAtoms) = Vy(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vy(nAtoms + 1:nAtoms + numAtoms)) + VY0;

nAtoms = nAtoms + numAtoms;

end
