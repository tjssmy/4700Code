global LAtoms WAtoms

%define globals to set number of atoms in rectangle
LAtoms = 10;        
WAtoms = 10;

doPlot = 1;
dt = 15e-15;
TStop = 1000 * dt;
InitDist = 0.0;
Method = 'VE'; % VE -- verlot; FD -- Forward Difference

Mass0 = 14 * C.am; % Silicon
Mass1 = 5 * C.am; % Argon

AtomSpacing = 0.5430710e-9;
LJSigma = AtomSpacing / 2^(1/6);
LJEpsilon = 1e-21;

PhiCutoff = 3 * AtomSpacing * 1.1;      %cut-off set to more than three times the atomic spacing - these atoms will not interact with one another

T = 30;                                 %temperature set to 30 degrees

AddRectAtomicArray(LAtoms, WAtoms, 0, 0, 0, 0, 0, T, 0);

%sets up the plot sizing
Size = 8 * AtomSpacing;
Limits = [-Size +Size -Size +Size]; % square is good
PlDelt = 5*dt;
MarkerSize = 10;
PlotFile = 'Block.gif';
doPlotImage = 1;
PlotSize = [100, 100, 1049, 1049];

ScaleV = .2e-11;
ScaleF = 20;
