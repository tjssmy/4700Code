doPlot = 1;
dt = 5e-15;
TStop = 3000 * dt;
InitDist = 0.0;
Method = 'VE'; % VE -- verlot; FD -- Forward Difference

Mass0 = 14 * C.am; % Silicon
Mass1 = 100 * C.am; % Argon

AtomSpacing = 0.5430710e-9;
LJSigma = AtomSpacing / 2^(1 / 6);
LJEpsilon = 1e-21;

PhiCutoff = 3 * AtomSpacing * 1.1;

T = 30;

LAtoms=10;
WAtoms=10;
X0=0;
Y0=0;
VX0=0; 
VY0=0;
InitDist=0; 
Temp=T;
Type=0;
AddRectAtomicArray(LAtoms, WAtoms, X0, Y0, VX0, VY0, InitDist, Temp, Type);

% vy0 = -sqrt(0.02*Ep/Mass1);
% AddRectAtomicArray(4,4,0,12*AtomSpacing,0,vy0,0,T,1);

Ep = 2;
numStream=5;
x0=0.1;
y0=10;
PartAng=-pi/2;
Type=1;
Ep=Ep*C.q_0;
Seper=5;
AddParticleStream(numStream, x0, y0, PartAng, Type, Ep, Seper);

Size = 10*AtomSpacing;
Limits = [-Size +Size -Size +Size]; % square is good
PlDelt = 5 * dt;

PlotFile = 'BlockSt.gif';
doPlotImage = 1;
PlotSize = [100, 100, 1049, 1049];

ScaleV = .02e-11;
ScaleF = 10;
