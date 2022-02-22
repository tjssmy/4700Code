Coupled = 1;
TwoCarriers = 0;
RC = 0;

nx = 201;
l = 1e-6;

x =linspace(0,l,nx);
dx = x(2)-x(1);
xm = x(1:nx-1) + 0.5*dx;

Nd = 1e16 * 1e6; % Const. 1/cm3 (100 cm/m)^3
NetDoping = ones(1,nx).*Nd; % doping

x0 = l/2;
nw = l/20;
npDisturbance = 1e16*1e6*exp(-((x-x0)/nw).^2);

LVbc = 1;
RVbc = 0;

TStop = 1200000*1e-18;
PlDelt = 10000*1e-18;

PlotSS = 1;

PlotYAxis = {[0 1] [-2e6 2e6] [-2e3 0]...
    [1e22 2e22] [0 1e22] [0 20e43]...
    [-5e34 5e34] [0 20e33] [-6e8 6e8] ...
    [0 20e7] [0 2e-2] [0 2e23]};

doPlotImage = 1;
PlotFile = 'Gau1Car1V.gif';
