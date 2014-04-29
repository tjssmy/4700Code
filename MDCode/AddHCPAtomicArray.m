function AddHCPAtomicArray(LAtoms,WAtoms,X0,Y0,VX0,VY0,RotAng,InitDist,Temp,Type)
global C
global x y AtomSpacing
global nAtoms  MinX MaxX MinY MaxY
global AtomType Vx Vy Mass0 Mass1

if Type == 0
    Mass = Mass0;
else
    Mass = Mass1;
end

L = ((LAtoms-1)+0.5)*AtomSpacing;
W = (WAtoms-1)*sqrt(3)/2*AtomSpacing;

p =1;
n = 4;
for i=1:LAtoms
    for j = 1:WAtoms
        y(end+1) = (sqrt(3)*(j-1)*AtomSpacing)/2.0;
        if rem(j,2) == 1
            x(end+1) = (i-1)*AtomSpacing;
        else
            x(end+1) = AtomSpacing/2.0 + (i-1)*AtomSpacing;
        end
        p = p + 1;
    end
end

x = x - L/2;
y = y - W/2;

numAtoms = p-1;

x(nAtoms+1:nAtoms+numAtoms) = x(nAtoms+1:nAtoms+numAtoms) + ...
    (rand(1,numAtoms)-0.5)*AtomSpacing*InitDist + X0;
y(nAtoms+1:nAtoms+numAtoms) = y(nAtoms+1:nAtoms+numAtoms) + ...
    (rand(1,numAtoms)-0.5)*AtomSpacing*InitDist + Y0;

AtomType(nAtoms+1:nAtoms+numAtoms) = Type;

if Temp == 0
    Vx(nAtoms+1:nAtoms+numAtoms) = 0;
    Vy(nAtoms+1:nAtoms+numAtoms) = 0;
else
    std0 = sqrt(C.kb*Temp/Mass);
    
    Vx(nAtoms+1:nAtoms+numAtoms) = std0*randn(1,numAtoms);
    Vy(nAtoms+1:nAtoms+numAtoms) = std0*randn(1,numAtoms);
end

Vx(nAtoms+1:nAtoms+numAtoms) = Vx(nAtoms+1:nAtoms+numAtoms) - ...
    mean(Vx(nAtoms+1:nAtoms+numAtoms)) + VX0;
Vy(nAtoms+1:nAtoms+numAtoms) = Vy(nAtoms+1:nAtoms+numAtoms) - ...
    mean(Vy(nAtoms+1:nAtoms+numAtoms)) + VY0;

nAtoms = nAtoms+numAtoms;

end

