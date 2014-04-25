function RectAtomicArray(LAtoms,WAtoms,Type)
global x y AtomSpacing 
global nAtoms  MinX MaxX MinY MaxY 
global AtomType

L = (LAtoms-1)*AtomSpacing;
W = (WAtoms-1)*AtomSpacing;

numAtoms = LAtoms*WAtoms;

xp(1,:) = linspace(0,L,LAtoms);
yp(1,:) = linspace(0,W,WAtoms);

x(nAtoms+1:nAtoms+LAtoms) = xp;
y(nAtoms+1:nAtoms+LAtoms) = yp(1);

for i = 1:WAtoms-1
    x(i*LAtoms+1:(i+1)*LAtoms) = xp;
    y(i*LAtoms+1:(i+1)*LAtoms) = yp(i+1);
end

x = x + (rand(1,nAtoms)-0.5)*AtomSpacing*InitDist;
y = y + (rand(1,nAtoms)-0.5)*AtomSpacing*InitDist;
AtomType(nAtoms+1:nAtoms+numAtoms) = Type;

MaxX = max(max(x)*1.1,MaxX);
MinX = max(mam(x)*1.1,MaxX);
nAtoms = nAtoms+numAtoms;

end

