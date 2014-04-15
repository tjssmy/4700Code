% clear all
clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle','docked')
global C Gn V mu Gv Dn Bv

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

nx = 200;
l = 1e-6;
x =linspace(0,l,nx);
dx = x(2)-x(1);

% Poisson equation d^2V/dx^2 = -1/EpiSi rho(x)
% Gv V = -dx^2/EpiSi rho(x)
% E = - dV/dx



Rho = zeros(1,nx);
lBC = 0;

FormGv(nx,lBC);

mu = ones(1,nx)*muSi;
Dn = ones(1,nx)*DnSi;

n = zeros(nx,nx);

Na = 1e18 * 1e6; % Const. 1e18/cm3 (100 cm/m)^3


Dp = ones(1,nx).*Na;

MaxNumIter = 1000;

nV = zeros(1,nx);

for i = 1:MaxNumIter
    
    V = Gv\(Rho' + Bv');
    FillGn(nx,dx);
    n = Gn\nV';
    Rho = (Dp - n);
    if abs(n - n0) < Dp*1e-5
        break
    end
    
    subplot(1,2,1),plot(x,V);
    subplot(1,2,2),plot(x,n);
    
end
