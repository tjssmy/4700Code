% clear all
clearvars
clearvars -GLOBAL
close all
format shorte
set(0,'DefaultFigureWindowStyle','docked')
global C Gn V mu Gv Dn Bv Em DnM muM A B np x xm n divFn
global Rho

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb*2*pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light

EpiSi = C.eps_0*11.68;
muSi = 450*1e-4; % cm2 V-1s-1 * 1/(100 cm/m)^2
DnSi = muSi*C.kb*300/C.q_0; % D = kt/q mu
niSi = 1e10*1e6; % 1/cm^3 * (100 cm/m)^3 intrinsic concentration

Coupled = 1;

nx = 201;
l = 1e-6;
x =linspace(0,l,nx);

dx = x(2)-x(1);
xm = x(1:nx-1) + 0.5*dx;

% Poisson equation d^2V/dx^2 = -1/EpiSi rho(x)
% Gv V = -dx^2/EpiSi rho(x)
% E = - dV/dx


FormGv(nx,0);
nV = zeros(1,nx);
[L,U] = lu(Gv);

mu = ones(1,nx)*muSi;
Dn = ones(1,nx)*DnSi;

muM(1:nx-1) = (mu(1:nx-1) + mu(2:nx))/2;
DnM(1:nx-1) = (Dn(1:nx-1) + Dn(2:nx))/2;
  
n = zeros(1,nx);

Nd = 1e16 * 1e6; % Const. 1/cm3 (100 cm/m)^3
Dp = ones(1,nx).*Nd; % doping
% n0 = (Dp + sqrt(Dp.^2 + 4* niSi*niSi))/2;

% Dp = zeros(1,nx);
x0 = l/2;
nw = l/20;
n0 = Dp + 1e16*1e6*exp(-((x-x0)/nw).^2);
divFn = zeros(1,nx);
dtMax = dx^2/2/max(Dn);

Rho = zeros(1,nx);
if (Coupled)
        Rho = C.q_0*(Dp - n0); % update Rho
        Rho(1) = 0;
end

Rho(1) = 0;
LVbc = 1;
Bv(1) = LVbc;

V = U\(L\(-dx^2/EpiSi*Rho' + Bv'));
Em(1:nx-1) = -(V(2:nx) - V(1:nx-1))/dx;
MaxEm = max(abs(Em));
Maxn = max(n0);
Ld = sqrt(EpiSi/(C.q_0*Maxn));
dxMax = Ld/5;

if MaxEm > 0
    dt = min([2*dx/MaxEm dtMax]);
else
    dt = dtMax
end

t = 0;
n = n0;
np = n0;

PlotVals(nx,dx,'on',[]);

TStop = 30000*dt;
PlDelt = 1000*dt;
Plt0 = PlDelt;
while t < TStop
    fprintf('time: %g\n',t);
    
    MaxEm = max(abs(Em));
    
    if MaxEm > 0
        dt = min([2*dx/MaxEm dtMax]);
    else
        dt = dtMax;
    end
    
    V = U\(L\(-dx^2/EpiSi*Rho' + Bv'));
    Em(1:nx-1) = -(V(2:nx) - V(1:nx-1))/dx;
        
    DivFn(nx,dx);
    n = np + dt*divFn;
    Maxn = max(n);
    
    Ld = sqrt(EpiSi/(C.q_0*Maxn));
    dxMax = Ld/5;

    if (Coupled)
        Rho = C.q_0*(Dp - n); % update Rho
        Rho(1) = 0;
    end
  
   
    np = n;
    t = t + dt;
    if t > Plt0
        PlotVals(nx,dx,'on',[]);
        Plt0 = Plt0 + PlDelt;
        pause(0.0001)
    end
    
  
end

figure
PlotVals(nx,dx,'off',[]);
