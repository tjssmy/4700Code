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
global MinX MaxX MinY MaxY PhiTot KETot
global nAtoms0 nAtoms1 T T0 T1 MarkerSize 
global doPlotImage PlotCount map im PlotSize

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
MarkerSize = 12;
Limits = [];
doPlot = 1;
doPlotImage = 0;
PlotCount  = 0;
PlotFile = 'image.gif';
PlotSize = [100, 100, 1049, 895];
% InitVStream
% InitHCP
% InitHCPBlob
InitVStreamHCP

MaxX = max(x)*1.5;
MinX = min(x)*1.5;

MaxY = max(y)*1.5;
MinY = min(y)*1.5;

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

V2 = Vx.*Vx + Vy.*Vy;
% KEc = 1/2*Mass*mean(V2);
% Tc = KEc/C.kb;

t = 0;
c = 1;
time(c) = 0;

PhiTot(c) = sum(Phi)/2;

V2_0 = (Vx(Pty0in).*Vx(Pty0in)+Vy(Pty0in).*Vy(Pty0in));
if nAtoms1
    V2_1 = (Vx(Pty1in).*Vx(Pty1in)+Vy(Pty1in).*Vy(Pty1in));
else 
    V2_1 = 0;
end

KE0 = mean(V2_0)*Mass0*0.5;
KE1 = mean(V2_1)*Mass1*0.5;

KETot(c) = (KE0*nAtoms0 + KE1*nAtoms1);
T(c) = KETot(c)/nAtoms/C.kb;
T0(c) = KE0/C.kb;
T1(c) = KE1/C.kb;

PlotVars(c,Limits);

xp = x - dt*Vx;
xpp = x - 2*dt*Vx;
yp = y - dt*Vy;
ypp = y - 2*dt*Vy;

Plt0 = PlDelt;

while t < TStop
    
    %     F = ma
    %     F = m dv/dt
    
    GetForces(PhiCutoff,LJEpsilon,LJSigma);
    
    % Forward difference
    if Method == 'FD'
        %     dv = F/m dt
        %     x = Vx*dt + F/m (dt)^2/2
        
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
    
    
    PhiTot(c) = sum(Phi)/2;
    V2_0 = (Vx(Pty0in).*Vx(Pty0in)+Vy(Pty0in).*Vy(Pty0in));
    if nAtoms1
        V2_1 = (Vx(Pty1in).*Vx(Pty1in)+Vy(Pty1in).*Vy(Pty1in));
    else 
        V2_1 = 0;
    end
    
    KE0 = mean(V2_0)*Mass0*0.5;
    KE1 = mean(V2_1)*Mass1*0.5;
    
    KETot(c) = (KE0*nAtoms0 + KE1*nAtoms1);
    T(c) = KETot(c)/nAtoms/C.kb;
    T0(c) = KE0/C.kb;
    T1(c) = KE1/C.kb;
    
    if t > Plt0
        fprintf('time: %g (%5.2g %%)\n',t,t/TStop*100);
        
        PlotVars(c,Limits);
        
        Plt0 = Plt0 + PlDelt;
        pause(0.00001)
    end
    
end


if doPlotImage
    imwrite(im,map,PlotFile,'DelayTime',0.2,'LoopCount',inf);
end


