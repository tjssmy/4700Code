% clear all
clearvars
clearvars -GLOBAL
close all
% set(0,'DefaultFigureWindowStyle','docked')
global C fc

addpath ../geom2d/geom2d

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light

nx = 100; %nx and ny are the bounds for the graphs made
ny = 50; 

Acond = 1;
Bcond = 10;

SimType = 'e';

Max = 5;
ncircs = 20;
doPlot = 1; % set to 1 if you would like to draw the plot

if SimType == 'c'
    n = 20;
    nSims = 1;
    Res = zeros(n, nSims);

    for k = 1:n
        Max(k) = k;
        fc = k;
        for i = 1:nSims
            % V = 1 --> R = 1/I
            Res(k,i) = 1/GetCurrents(ncircs,Max(k),nx,ny,...
                Acond,Bcond,doPlot,SimType);
        end
    end
else %if SimType == 'e'

    n = 20;
    nSims = 30;
    Res = zeros(n, nSims);

    for k = 1:n
        Max(k) = k * 3;
        fc = k;
        for i = 1:nSims
            % V = 1 --> R = 1/I
            Res(k,i) = 1 / GetCurrents(Max(k), 10, nx, ny,...
                Acond, Bcond, doPlot, SimType);
        end
    end

end

if doPlot
    imwrite(im, map, 'imagefile.gif', 'DelayTime', 0.2, 'LoopCount', inf);
    figure
end

subplot(4, 1, 1), bar(Res(1, :))
xl = sprintf('Simulation Number (Num = %i)', Max(1));
xlabel(xl);
ylabel('Resistance');

subplot(4, 1, 2), bar(Res(n, :))
xl = sprintf('Simulation Number (Num = %i)', Max(n));
xlabel(xl);
ylabel('Resistance');

AveCurr = mean(Res');
StdCurr = std(Res');

subplot(4, 1, 3), plot(Max, AveCurr);
hold on
subplot(4, 1, 3), plot(Max, AveCurr - StdCurr, 'r');
subplot(4, 1, 3), plot(Max, AveCurr + StdCurr, 'r');
xlabel('Num');
ylabel('Resistance');
hold off

subplot(4, 1, 4), plot(Max, StdCurr./ AveCurr * 100, 'r');
xlabel('Num');
ylabel('Variance (% of mean)');
