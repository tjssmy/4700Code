% clear all
clearvars
clearvars -GLOBAL
close all
format shorte

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

doPlot = 1;
dt = 5e-15;
TStop = 10000*dt;
InitDist = 0.0;
Method = 'VE'; % VE -- verlot; FD -- Forward Difference

AddParticle = 4;
Ep = 10*C.q_0;

Mass = 14*C.am; % Silicon
AtomSpacing = 0.5430710e-9;
LJSigma = AtomSpacing/2^(1/6);
LJEpsilon = 1e-21;

LAtoms = 10;
WAtoms = 10;

L = (LAtoms-1)*AtomSpacing*1.0;
W = (WAtoms-1)*AtomSpacing*1.0;
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

if AddParticle
    PartAng = -pi/8;
    x0 = 1.57876*AtomSpacing;
    y0 = W+2.0*AtomSpacing;
    for p = 0:AddParticle-1
        nAtoms = nAtoms+1;
        x(nAtoms) = x0 - 3*p*AtomSpacing*cos(PartAng);
        y(nAtoms) = y0 - 3*p*AtomSpacing*sin(PartAng);
    end
end

Fx = zeros(1,nAtoms);
Fy = zeros(1,nAtoms);
Phi = zeros(1,nAtoms);

GetForces(PhiCutoff,LJEpsilon,LJSigma);
Phi0 = sum(Phi)/nAtoms/2;
KE = Phi0 + LJEpsilon;
T = KE/C.kb;
T = 30;
if T == 0
    Vx(1:nAtoms) = 0;
    Vy(1:nAtoms) = 0;
else
    std = sqrt(C.kb*T/Mass);
    Vx = std*randn(1,nAtoms);
    Vy = std*randn(1,nAtoms);
end

Vx = Vx - mean(Vx);
Vy = Vy - mean(Vy);

V2 = Vx.*Vx + Vy.*Vy;
KEc = 1/2*Mass*mean(V2);
Tc = KEc/C.kb;

if AddParticle
    V = sqrt(2*Ep/Mass);
    for p = 1:AddParticle
        Vx(nAtoms-AddParticle+p) = V*cos(PartAng);
        Vy(nAtoms-AddParticle+p) = V*sin(PartAng);
    end
end

V2 = Vx.*Vx + Vy.*Vy;
MaxV = max(sqrt(V2));
ScaleV = sqrt(L*L + W*W)/MaxV*0.5;

subplot(2,2,1),plot(x,y,'o','markers',20);
hold on
subplot(2,2,1),quiver(x,y,Fx,Fy,0);
hold off

axis([-2*AtomSpacing,(LAtoms+2)*AtomSpacing,...
    -2*AtomSpacing,(WAtoms+2)*AtomSpacing]);

t = 0;
c = 1;

PhiTot(c) = sum(Phi)/2;
time(c) = 0;
KETot(c) = KE*nAtoms;

subplot(2,2,3),plot(x,y,'o','markers',20);
hold on
subplot(2,2,3),quiver(x,y,Vx*ScaleV,Vy*ScaleV,0,'r');
hold off
axis([-2*AtomSpacing,(LAtoms+2)*AtomSpacing,...
    -2*AtomSpacing,(WAtoms+2)*AtomSpacing]);

subplot(2,2,2),plot(time,KETot);
hold on
subplot(2,2,2),plot(time,PhiTot+LJEpsilon*nAtoms,'g');
subplot(2,2,2),plot(time,PhiTot+KETot+LJEpsilon*nAtoms,'k');
%     axis([0 TStop 0e-22 3e-23]);

AE = 1/2*Mass*V2 + Phi-Phi0;
if AddParticle, AE(end) = 0;end
% AE = Phi;
subplot(2,2,4),scatter3(x,y,AE,ones(1,nAtoms)*300,AE);
axis([-2*AtomSpacing,(LAtoms+2)*AtomSpacing,...
    -2*AtomSpacing,(WAtoms+2)*AtomSpacing]);
view(2);

xp = x - dt*Vx;
xpp = x - 2*dt*Vx;
yp = y - dt*Vy;
ypp = y - 2*dt*Vy;

while t < TStop
    
    %     F = ma
    %     F = m dv/dt
    %     dv = F/m dt
    %     x = Vx*dt + F/m (dt)^2/2
    time(c) = t;
    GetForces(PhiCutoff,LJEpsilon,LJSigma);
    
    % Forward difference
    if Method == 'FD'
        dvx = Fx*dt/Mass;
        Vx = Vx + dvx;
        dx = Vx*dt + Fx*dt^2/2/Mass;
        
        dvy = Fy*dt/Mass;
        Vy = Vy + dvy;
        dy = Vy*dt + Fy*dt^2/2/Mass;
        
        x = xp + dx;
        y = yp + dy;
    elseif Method == 'VE'
        
        x = -xpp + 2*xp + dt^2/Mass*Fx;
        y = -ypp + 2*yp + dt^2/Mass*Fy;
        
        Vx = (x - xpp)/(2*dt);
        Vy = (y - ypp)/(2*dt);
    end
    
    V2 = Vx.*Vx + Vy.*Vy;
    KE = 1/2*Mass*mean(V2);
    KETot(c) = KE*nAtoms;
    
    T(c) = KE/C.kb;
    
    PhiTot(c) = sum(Phi)/2; % every pair contributes twice
    
    subplot(2,2,1),plot(x,y,'o','markers',20);
    hold on
    subplot(2,2,1),quiver(x,y,Fx,Fy,0);
    hold off
    axis([-2*AtomSpacing,(LAtoms+2)*AtomSpacing,...
        -2*AtomSpacing,(WAtoms+2)*AtomSpacing]);
    
    
    subplot(2,2,2),plot(time,KETot);
    hold on
    subplot(2,2,2),plot(time,PhiTot-Phi0*nAtoms,'g');
    subplot(2,2,2),plot(time,PhiTot-Phi0*nAtoms+KETot,'k');
    %     axis([0 TStop -0e-19 0.5e-19]);
    
    subplot(2,2,3),plot(x,y,'o','markers',20);
    hold on
    subplot(2,2,3),quiver(x,y,Vx*ScaleV,Vy*ScaleV,0,'r');
    hold off
    axis([-2*AtomSpacing,(LAtoms+2)*AtomSpacing,...
        -2*AtomSpacing,(WAtoms+2)*AtomSpacing]);
    
    AE = 1/2*Mass*V2 + Phi-Phi0;
    if AddParticle, AE(end) = 0;end
    % AE = Phi;
    subplot(2,2,4),scatter3(x,y,AE,ones(1,nAtoms)*300,AE);
    axis([-2*AtomSpacing,(LAtoms+2)*AtomSpacing,...
        -2*AtomSpacing,(WAtoms+2)*AtomSpacing]);
    view(2);
    
    pause(0.00001)
    
    xpp = xp;
    ypp = yp;
    
    xp = x;
    yp = y;
    
    c = c + 1;
    t  = t + dt;
end
