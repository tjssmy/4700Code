% clear all
clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle','docked')
global C im fig map fc

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

nTime = 100;
nTraj = 10000;
nSims = 100;

InitalAngle = 55*pi/180;
PlotTraj = 1;
MaxC = 10000;



V0 = 1000;
g = 1
c = 2;
dt = 1;

% Wind = @UniformRandWind;
% WindParas = [10];

% Wind = @NormalRandWind;
% WindParas = [1];

Wind = @ComplexRandWind;
WindParas = [5];


xl =[];

for n = 1: nSims
    x(1,:) = zeros(1,nTraj);
    y(1,:) = zeros(1,nTraj);
    
    Vx(1:nTraj) = V0*cos(InitalAngle);
    Vy(1:nTraj) = V0*sin(InitalAngle);

    for c=2:MaxC
        
        dvx = Wind(nTraj,WindParas);
        Vx = Vx + dvx;
        dx = Vx*dt;
        
        dvy = -g*dt;
        Vy = Vy + dvy;
        dy = Vy*dt + g*dt^2/2;
        
        x(c,:) = x(c-1,:) + dx;
        y(c,:) = y(c-1,:) + dy;
        if max(y(c)) < 0
            break
        end
    end
    
    subplot(2,1,1),plot(x,y);
    axis([0 13e5 0 4e5]);
    hold on
    
    xl = [xl x(end,:)];
    subplot(2,1,2),hist(xl,100)
    axis([0 13e5 0 50]);
    pause(0.01);
end