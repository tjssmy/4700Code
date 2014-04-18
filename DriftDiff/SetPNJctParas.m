
Coupled = 1;
TwoCarriers = 1;
RC = 1;

nx = 501;
l = 1e-6;

x =linspace(0,l,nx);
dx = x(2)-x(1);
xm = x(1:nx-1) + 0.5*dx;

ni = x < l/2;
pi = x >= l/2;

Nd = 2e16 * 1e6; % Const. 1/cm3 (100 cm/m)^3
Na = 1e16 * 1e6;
NetDoping(ni) = Nd;
NetDoping(pi) = -Na;

x0 = l/2;
nw = l/20;
npDisturbance = 0e16*1e6*exp(-((x-x0)/nw).^2);

LVbc = 'fl';
RVbc = 0;

TStop = 42000000*1e-18;
PlDelt = 200000*1e-18;

PlotYAxis = {[0 0.3] [0e5 15e5] [-20e2 30e2]...
    [0e21 2.5e22] [0 1.1e22] [0 20e43]...
    [-5e33 5e33] [-5e33 5e33] [-0e8 15e8] ...
    [1e-3 5e8] [0e-3 30e-3] [0 2.5e22]};