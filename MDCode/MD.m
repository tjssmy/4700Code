% clear all
clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle','docked')
global C x y nAtoms Fx Fy Phi

addpath ../geom2d/geom2d

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb*2*pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665; %metres (32.1740 ft) per s²
C.am = 1.66053892e-27;

nTime = 100;

doPlot = 1;
dt = 5e-15;
TStop = 10000*dt;
InitDist = 0.00;

Mass = 14*C.am; % Silicon
AtomSpacing = 0.5430710e-9;
LJSigma = AtomSpacing/2^(1/6);
LJEpsilon = 1e-21;

LAtoms = 1;
WAtoms = 2;

L = LAtoms*AtomSpacing*0.90;
W = WAtoms*AtomSpacing*0.90;
PhiCutoff = 3*AtomSpacing*1.1;

nAtoms = LAtoms*WAtoms;

xp(1,:) = linspace(0,L,LAtoms);
yp(1,:) = linspace(0,W,WAtoms);

x = xp;
y(1:LAtoms) = yp(1);

for i = 1:WAtoms-1
    x(i*LAtoms+1:(i+1)*LAtoms) = xp;
    y(i*LAtoms+1:(i+1)*LAtoms) = yp(i+1);
end

x = x + (rand(1,nAtoms)-0.5)*AtomSpacing*InitDist;
y = y + (rand(1,nAtoms)-0.5)*AtomSpacing*InitDist;

Fx = zeros(1,nAtoms);
Fy = zeros(1,nAtoms);
Phi = zeros(1,nAtoms);

Vx(1:nAtoms) = 0;
Vy(1:nAtoms) = 0;

subplot(2,2,1),plot(x,y,'o','markers',20);
hold on
subplot(2,2,1),quiver(x,y,Fx,Fy);
hold off

axis([-2*AtomSpacing,(LAtoms+2)*AtomSpacing,...
    -2*AtomSpacing,(WAtoms+2)*AtomSpacing]);

t = 0;
c = 1; 
while t < TStop
    
%     F = ma
%     F = m dv/dt
%     dv = F/m dt
%     x = Vx*dt + F/m (dt)^2/2
    time(c) = t;
    GetForces(PhiCutoff,LJEpsilon,LJSigma);
    
    dvx = Fx*dt/Mass;
    Vx = Vx + dvx;
    dx = Vx*dt + Fx*dt^2/2/Mass;
    
    dvy = Fy*dt/Mass;
    Vy = Vy + dvy;
    dy = Vy*dt + Fy*dt^2/2/Mass;
    
    x = x + dx;
    y = y + dy;
    
    V2 = Vx.*Vx + Vy.*Vy;
    KE(c) = 1/2*Mass*sum(V2);
    PhiTot(c) = sum(Phi);
    
    subplot(2,2,1),plot(x,y,'o','markers',20);
    hold on
    subplot(2,2,1),quiver(x,y,Fx,Fy);
    hold off
    axis([-2*AtomSpacing,(LAtoms+2)*AtomSpacing,...
        -2*AtomSpacing,(WAtoms+2)*AtomSpacing]);
    
    subplot(2,2,3),plot(x,y,'o','markers',20);
    hold on
    subplot(2,2,3),quiver(x,y,Vx,Vy,'r');
    hold off
    axis([-2*AtomSpacing,(LAtoms+2)*AtomSpacing,...
        -2*AtomSpacing,(WAtoms+2)*AtomSpacing]);

    subplot(2,2,2),plot(time,KE);
    hold on
%     axis([0 TStop -20e-19 3e-19]);
    subplot(2,2,2),plot(time,PhiTot,'g');
    subplot(2,2,2),plot(time,PhiTot+KE,'k');
    pause(0.001)
    
    c = c + 1;
    t  = t + dt;
end
