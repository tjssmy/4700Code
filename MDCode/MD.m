% clear all
clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle','docked')
global C x y nAtoms Fx Fx

addpath ../geom2d/geom2d

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb*2*pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665 %metres (32.1740 ft) per s²
C.am = 1.66053892e10-27;

nTime = 100;

doPlot = 1;
TStop = 100e-15;
InitDist = 0.1;
Mass = 14*C.am; % Silicon
AtomSpacing = 0.5430710e-9;

LAtoms = 10;
WAtoms = 8;

L = LAtoms*AtomSpacing;
W = WAtoms*AtomSpacing;
PhiCutoff = 3*AtomSpacing;

nAtoms = LAtoms*WAtoms;

dt = 1e-15;

xp(1,:) = linspace(0,L,LAtoms);
yp(1,:) = linspace(0,W,WAtoms);

x = xp;
y(1:LAtoms) = yp(1);

for i = 1:WAtoms-1
    x(i*LAtoms+1:(i+1)*LAtoms) = xp;
    y(i*LAtoms+1:(i+1)*LAtoms) = yp(i+1);
end

x = x + rand(1,nAtoms)*AtomSpacing*InitDist;
y = y + rand(1,nAtoms)*AtomSpacing*InitDist;

Fx = zeros(1,nAtoms);
Fy = zeros(1,nAtoms);

Vx = Fx*dt;
Vy = Fy*dt;

plot(x,y,'o','markers',48);
hold
quiver(x,y,Vx,Vy);

axis([-2*AtomSpacing,(LAtoms+2)*AtomSpacing,...
    -2*AtomSpacing,(WAtoms+2)*AtomSpacing]);

t = 0;

while t < TStop
    
%     dvx = *dt;
    Vx = Vx + dvx;
    dx = Vx*dt;
    
    dvy = -g*dt;
    Vy = Vy + dvy;
    dy = Vy*dt + g*dt^2/2;
    
    x(c,:) = x(c-1,:) + dx;
    y(c,:) = y(c-1,:) + dy;
    
    t  = t + dt;
end
