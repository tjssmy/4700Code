function AddRectAtomicArray(LAtoms, WAtoms, X0, Y0, VX0, VY0, InitDist, Temp, Type)
global C
global x y AtomSpacing
global nAtoms
global AtomType Vx Vy Mass0 Mass1

if Type == 0            %sets mas based on the type of atom set in MD.m (global variable)
    Mass = Mass0;
else
    Mass = Mass1;
end

%sets the lenght and with based on the number of atoms provided when
%calling this function
L = (LAtoms - 1) * AtomSpacing;
W = (WAtoms - 1) * AtomSpacing;

numAtoms = LAtoms * WAtoms;

xp(1, :) = linspace(0, L, LAtoms);
yp(1, :) = linspace(0, W, WAtoms);

x(nAtoms + 1:nAtoms+LAtoms) = xp-L/2;
y(nAtoms + 1:nAtoms+LAtoms) = yp(1)-W/2;

for i = 1:WAtoms-1
    x(nAtoms + i * LAtoms + 1:nAtoms + (i + 1) * LAtoms) = xp - L / 2;
    y(nAtoms + i * LAtoms + 1:nAtoms + (i + 1) * LAtoms) = yp(i + 1) - W / 2;
end

x(nAtoms + 1:nAtoms + numAtoms) = x(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + X0;
y(nAtoms + 1:nAtoms + numAtoms) = y(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + Y0;

AtomType(nAtoms + 1:nAtoms + numAtoms) = Type;

if Temp == 0
    Vx(nAtoms + 1:nAtoms + numAtoms) = 0;
    Vy(nAtoms + 1:nAtoms + numAtoms) = 0;
else
    std0 = sqrt(C.kb * Temp / Mass);

    Vx(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
    Vy(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
end

Vx(nAtoms + 1:nAtoms + numAtoms) = Vx(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vx(nAtoms + 1:nAtoms + numAtoms)) + VX0;
Vy(nAtoms + 1:nAtoms + numAtoms) = Vy(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vy(nAtoms + 1:nAtoms + numAtoms)) + VY0;

nAtoms = nAtoms + numAtoms;

end
