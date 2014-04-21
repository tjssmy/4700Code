
Coupled = 1;
TwoCarriers = 1;
RC = 1;

nx = 201;
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

JBC = 1;

RVbc = 0;

TStop = 8000000*1e-18;
PlDelt = 400000*1e-18;

Phi =  C.Vt *log(Na*Nd/(niSi*niSi));
W  = sqrt(2*EpiSi*(Nd+Nd)*(Phi)/(C.q_0*Nd*Na));
Wn = W*Na/(Nd+Na);
Wp = (W - Wn);

LVbc = Phi;


PlotYAxis = {[0 Phi+0.1] [0e5 40e5] [-20e2 40e2]...
    [0e21 2.5e22] [0 1.1e22] [0 20e43]...
    [-5e33 5e33] [-5e33 5e33] [-0e8 3e8] ...
    [1e-3 1e8] [0e-3 30e-3] [0 2.5e22]};

PlotFile = 'PNJuncEq.gif';
doPlotImage = 1


SecondSim = 1;
LVbc2 = 'fl';
TStop2 = TStop +  8000000*1e-18;

fprintf('Phi: %g W: %g \n',Phi,W)

