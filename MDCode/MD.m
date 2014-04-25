% clear all
clearvars
clearvars -GLOBAL
close all
format shorte

set(0,'DefaultFigureWindowStyle','docked')
global C 
global Vx Vy x y Fx Fy AtomSpacing 
global Phi nAtoms time Mass0 Mass1 Pty0in Pty1in
global LJEpsilon LJSigma Phi0 AtomType
global MinX MaxX MinY MaxY

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

MaxX = 0;
MinX = 0;
MaxY = 0;
MinY = 0;
nAtoms = 0;

doPlot = 1;
dt = 5e-15;
TStop = 10000*dt;
InitDist = 0.0;
Method = 'VE'; % VE -- verlot; FD -- Forward Difference

Mass0 = 14*C.am; % Silicon
Mass1 = 250*C.am; % Argon

AddParticle = 1;
Ep = 0.5*C.q_0;
AddParticleType = 1;
AddParticleMass = Mass1;

AtomSpacing = 0.5430710e-9;
LJSigma = AtomSpacing/2^(1/6);
LJEpsilon = 1e-21;

PhiCutoff = 3*AtomSpacing*1.1;

RectAtomicArray(LAtoms,WAtoms,0)


AtomType(1:nAtoms) = 0;

if AddParticle
    PartAng = -pi/2;
    x0 = 5.57876*AtomSpacing;
    y0 = W+2.0*AtomSpacing;
    for p = 0:AddParticle-1
        nAtoms = nAtoms+1;
        x(nAtoms) = x0 - 3*p*AtomSpacing*cos(PartAng);
        y(nAtoms) = y0 - 3*p*AtomSpacing*sin(PartAng);
        AtomType(nAtoms) = AddParticleType;
    end
end

Fx = zeros(1,nAtoms);
Fy = zeros(1,nAtoms);
Phi = zeros(1,nAtoms);
dx = zeros(1,nAtoms);
dy = zeros(1,nAtoms);
dvx = zeros(1,nAtoms);
dvy = zeros(1,nAtoms);

Pty0in = AtomType == 0;
Pty1in = AtomType == 1;

nAtoms0 = sum(Pty0in);
nAtoms1 = sum(Pty1in);

GetForces(PhiCutoff,LJEpsilon,LJSigma);
Phi0 = sum(Phi)/nAtoms/2;
KE = Phi0 + LJEpsilon;
T = KE/C.kb;
T = 30;
if T == 0
    Vx(1:nAtoms) = 0;
    Vy(1:nAtoms) = 0;
else
    std0 = sqrt(C.kb*T/Mass0);
    std1 = sqrt(C.kb*T/Mass1);
    
    Vx(Pty0in) = std0*randn(1,nAtoms0);
    Vy(Pty0in) = std0*randn(1,nAtoms0);
    
    Vx(Pty1in) = std1*randn(1,nAtoms1);
    Vy(Pty1in) = std1*randn(1,nAtoms1);
    
end

Vx = Vx - mean(Vx);
Vy = Vy - mean(Vy);

V2 = Vx.*Vx + Vy.*Vy;
% KEc = 1/2*Mass*mean(V2);
% Tc = KEc/C.kb;

if AddParticle
    V = sqrt(2*Ep/AddParticleMass);
    for p = 1:AddParticle
        Vx(nAtoms-AddParticle+p) = V*cos(PartAng);
        Vy(nAtoms-AddParticle+p) = V*sin(PartAng);
    end
end

t = 0;
c = 1;
time(c) = 0;

PlotVars(c,AddParticle);

xp = x - dt*Vx;
xpp = x - 2*dt*Vx;
yp = y - dt*Vy;
ypp = y - 2*dt*Vy;

while t < TStop
    
    %     F = ma
    %     F = m dv/dt
    %     dv = F/m dt
    %     x = Vx*dt + F/m (dt)^2/2
   
    GetForces(PhiCutoff,LJEpsilon,LJSigma);
    
    % Forward difference
    if Method == 'FD'
        dvx(Pty0in) = Fx(Pty0in)*dt/Mass0;
        dvx(Pty1in) = Fx(Pty1in)*dt/Mass1;
        
        Vx = Vx + dvx;
        dx(Pty0in) = Vx(Pty0in)*dt + Fx(Pty0in)*dt^2/2/Mass0;
        dx(Pty1in) = Vx(Pty1in)*dt + Fx(Pty1in)*dt^2/2/Mass1;
        
        dvy(Pty0in) = Fy(Pty0in)*dt/Mass0;
        dvy(Pty1in) = Fy(Pty1in)*dt/Mass1;
        
        Vy = Vy + dvy;
        
        dy(Pty0in) = Vy(Pty0in)*dt + Fy(Pty0in)*dt^2/2/Mass0;
        dy(Pty1in) = Vy(Pty1in)*dt + Fy(Pty1in)*dt^2/2/Mass1;
        
        dy = Vy*dt + Fy*dt^2/2/Mass;
        
        x = xp + dx;
        y = yp + dy;
        
    elseif Method == 'VE'
        
        x(Pty0in) = -xpp(Pty0in) + 2*xp(Pty0in) + dt^2/Mass0*Fx(Pty0in);
        x(Pty1in) = -xpp(Pty1in) + 2*xp(Pty1in) + dt^2/Mass1*Fx(Pty1in);
        
        y(Pty0in) = -ypp(Pty0in) + 2*yp(Pty0in) + dt^2/Mass0*Fy(Pty0in);
        y(Pty1in) = -ypp(Pty1in) + 2*yp(Pty1in) + dt^2/Mass1*Fy(Pty1in);
        
        Vx = (x - xpp)/(2*dt);
        Vy = (y - ypp)/(2*dt);
    end
    
    xpp = xp;
    ypp = yp;
    
    xp = x;
    yp = y;
    
    c = c + 1;
    t  = t + dt;
    time(c) = t;
     
    PlotVars(c,AddParticle);
    pause(0.00001)
    
end
