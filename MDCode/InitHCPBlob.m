doPlot = 1;
dt = 15e-15;
TStop = 1000 * dt;
InitDist = 0.05;
Method = 'VE'; % VE -- verlot; FD -- Forward Difference

Mass0 = 14 * C.am; % Silicon
Mass1 = 4 * C.am; % Argon

AtomSpacing = 0.5430710e-9;
LJSigma = AtomSpacing / 2^(1 / 6);
LJEpsilon = 1e-21;

PhiCutoff = 3 * AtomSpacing * 1.1;

T = 30;

AddHCPAtomicBlob(10, 0, 0, 0, 0, 0, InitDist, T, 0);

Ep = 4 * C.q_0;
vy0 = -sqrt(0.02 * Ep / Mass1);
AddHCPAtomicBlob(5, -0, 20, 0, vy0, 0, InitDist, T, 1);

Size = 15 * AtomSpacing;
Limits = [-Size +Size -Size +1.5*Size]; % square is good
PlDelt = 5 * dt;
MarkerSize = 4;
PlotFile = 'HCPBlob.gif';
doPlotImage = 1;
PlotSize = [100, 100, 1049, 1549];

ScaleV = .2e-11;
ScaleF = 20;
