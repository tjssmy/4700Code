% clear all
% clearvars
close all
% set(0,'DefaultFigureWindowStyle','docked')
global C im fig map fc

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
doPlot = 1;


for k = 1:2
    MaxRad(k) = k*2
    fc = k;
    for i = 1:1
        % V = 1 --> R = 1/I
        Res(k,i)  = 1/GetCurrents(ncircs,MaxRad(k),nx,ny,Acond,Bcond,doPlot);
    end
end

if doPlot
    imwrite(im,map,'imagefile.gif','DelayTime',100,'LoopCount',inf);
end

subplot(4,1,1),bar(Res(1,:))
xl = sprintf('Simulation Number (MaxRad = %i)',MaxRad(1));
xlabel(xl);
ylabel('Resistance');

subplot(4,1,2),bar(Res(10,:))
xl = sprintf('Simulation Number (MaxRad = %i)',MaxRad(10));
xlabel(xl);
ylabel('Resistance');

% end

AveCurr = mean(Res');
StdCurr = std(Res');

subplot(4,1,3),plot(MaxRad,AveCurr);
hold on
subplot(4,1,3),plot(MaxRad,AveCurr-StdCurr,'r');
subplot(4,1,3),plot(MaxRad,AveCurr+StdCurr,'r');
xlabel('MaxRad');
ylabel('Resistance');
hold off

subplot(4,1,4),plot(MaxRad,StdCurr./AveCurr*100,'r');
xlabel('MaxRad');
ylabel('Variance (% of mean)');

