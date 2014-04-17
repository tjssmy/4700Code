% clear all
clearvars
clearvars -GLOBAL
close all
format shorte
set(0,'DefaultFigureWindowStyle','docked')
global C Gn V mu Gv Dn Bv Em DnM muM A B np

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb*2*pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light

EpiSi = C.eps_0*11.68
muSi = 450*1e-4 % cm2 V-1s-1 * 1/(100 cm/m)^2
DnSi = muSi*C.kb*300/C.q_0; % D = kt/q mu
niSi = 1e10*1e6; % 1/cm^3 * (100 cm/m)^3 intrinsic concentration

Coupled = 1;
SFG = 0;

nx = 201;
l = 1e-6;
x =linspace(0,l,nx);

dx = x(2)-x(1);
xm = x(1:nx-1) + 0.5*dx;

% Poisson equation d^2V/dx^2 = -1/EpiSi rho(x)
% Gv V = -dx^2/EpiSi rho(x)
% E = - dV/dx


Rho = zeros(1,nx);


FormGv(nx,0);
nV = zeros(1,nx);
[L,U,P,Q] = lu(Gv);

mu = ones(1,nx)*muSi;
Dn = ones(1,nx)*DnSi;

n = zeros(nx,nx);

Nd = 2e14 * 1e6; % Const. 1/cm3 (100 cm/m)^3

MaxNumIter = 1000;

Dp = ones(1,nx).*Nd; % doping
n0 = (Dp + sqrt(Dp.^2 + 4* niSi*niSi))/2;
nV(1) = n0(1);      % n0 at left end

lbc = -1e-3;

Bv(1) = lbc;

dtMax = dx^2/2/max(Dn);
dt = 

TStop = 1e-9;

Rho = C.q_0*(Dp - n0');
Rho(1) = 0;
V = U\(L\(-dx^2/EpiSi*Rho' + Bv'));
Em(1:nx-1) = -(V(2:nx) - V(1:nx-1))/dx;
MaxEm = max(abs(Em));
Maxn = max(n0);
Ld = sqrt(EpiSi/(C.q_0*Maxn));

t = 0;
n = n0;
np = n;

while t < TStop
    MaxEm = max(abs(Em));
    if MaxEM > 0
        dt = min([2*dx/MaxEM dtMax]);
    else
        dt = dtMax
    end
    

    V = U\(L\(-dx^2/EpiSi*Rho' + Bv'));
    
    divJn = DivJn(nx,dx);
    n = np - dt*divJn;
    
    if (Coupled)
        Rho = C.q_0*(Dp - n'); % update Rho
        Rho(1) = 0;
    end
    
    gradn = (n(2:nx)-n(1:nx-1))/dx;
    nM = (n(2:nx)+n(1:nx-1))/2;
    
    JnDiff = C.q_0*DnM.*gradn';
    JnDrift = C.q_0*muM.*nM'.*Em;
    
    subplot(2,2,1),plot(x,V);
    %     axis([0 l -0.05 0.01])
    hold on
    subplot(2,2,2),plot(xm,Em);
    %     axis([0 l -15e4 1e4])
    hold on
    subplot(2,2,3),plot(x,n);
    hold on
    %     axis([0 l 0 5e20])
    subplot(2,2,4),plot(xm,JnDiff,'r');
    hold on
    subplot(2,2,4),plot(xm,JnDrift,'b');
    %     axis([0 l -5e5 5e5])
    
    np = n;
    
    pause(0.0001)
end

