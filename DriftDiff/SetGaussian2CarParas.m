
Coupled = 1;
TwoCarriers = 1;
RC = 1;

nx = 2001;
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

TStop = 2200000*1e-18;
PlDelt = 20000*1e-18;

PlotYAxis = {[0 1] [0 2e6] [-1e4 1e4]...
    [0 2e22] [0 1e22] [0 20e43]...
    [-5e34 15e34] [-20e33 20e33] [-1e8 4e8] ...
    [-20e7 20e7] [0 20e-3] [0 2e23]};