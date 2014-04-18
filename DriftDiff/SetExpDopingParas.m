
Coupled = 1;
TwoCarriers = 1;
RC = 1;

nx = 201;
l = 1e-6;

x =linspace(0,l,nx);
dx = x(2)-x(1);
xm = x(1:nx-1) + 0.5*dx;

Nd = 2e16 * 1e6; % Const. 1/cm3 (100 cm/m)^3
a = 5e6;
NetDoping = Nd*exp(-a.*x); % doping

x0 = l/2;
nw = l/20;
npDisturbance = 0e16*1e6*exp(-((x-x0)/nw).^2);

LVbc = 'fl';
RVbc = 0;

TStop = 42000000*1e-18;
PlDelt = 200000*1e-18;

PlotYAxis = {[0 .1] [-2e5 2e5] [-1.5e2 1.5e2]...
    [0.5e21 2e22] [0 1e12] [0 20e32]...
    [-20e32 15e32] [-5e22 5e22] [-0e8 0.5e8] ...
    [-0e-3 1e-3] [0e-3 20e-3] [0 2e22]};