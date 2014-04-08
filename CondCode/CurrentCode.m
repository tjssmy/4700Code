% clear all
clearvars
close all
set(0,'DefaultFigureWindowStyle','docked')
global C

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb*2*pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light




nx = 100;
ny = 50;

Acond = 1;
Bcond = 10;

MaxRad = 5;
ncircs = 20;
doPlot = 0;

for k = 1:10
    MaxRad = k*2
    for i = 1:30
        Curr(k,i)  = GetCurrents(ncircs,MaxRad,nx,ny,Acond,Bcond,doPlot);
    end
end

subplot(2,1,1),bar(Curr(1,:))
subplot(2,1,2),bar(Curr(10,:))

% end

AveCurr = mean(Curr');
StdCurr = stdev(Curr');