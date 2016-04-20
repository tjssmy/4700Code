clear all
clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle','docked')

global C

addpath ../geom2d/geom2d

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665 %metres (32.1740 ft) per s²

nTime = 100;
nTraj = 5;
nSims = 100;

InitalAngle = 55 * pi / 180;
MaxC = 10000;


V0 = 1000;
g = 1;
c = 2;
dt = 1;

Wind = @UniformRandWind;
WindParas = [10];

% Wind = @NormalRandWind;
% WindParas = [5];

% Wind = @ComplexRandWind;
% WindParas = [.35];

m = 10;

% Added Energy to the trajectories using the equation E = Ek + Eg

E = 0.5 * m .* (V0 * sin(InitalAngle))^2;
initialEnergy = mean(E); % Initial Mean Energy of the trajectories;

Ek = zeros(1, nTraj);   %Set kinetic energy to zero
Eg = zeros(1, nTraj);   %Set gravitational potential energy to zero


for n = 1: nSims
    x(1, :) = zeros(1, nTraj);
    y(1, :) = zeros(1, nTraj);

    Vx(1:nTraj) = V0 * cos(InitalAngle);
    Vy(1:nTraj) = V0 * sin(InitalAngle);

    for c=2:MaxC

        dvx = Wind(nTraj,WindParas)*dt;
        Vx = Vx + dvx;
        dx = Vx * dt;

        dvy = -g * dt;
        Vy = Vy + dvy;
        dy = Vy * dt + g * dt^2 / 2;

        x(c,:) = x(c - 1,:) + dx;
        y(c,:) = y(c - 1,:) + dy;
        
        Eg = m .* g .* y(c,:);  % Calculate gravitational potential energy as:
                                % Eg = m*g*h where h = y
                                
        randLoss = (0.75 + 0.2 .* rand(1, nTraj)); % Random energy loss in the system 
                                                   % between 75-95%.
        
        Ek = E - randLoss .* Eg; % Calculate kinetic energy as:
                                   % Ek = E - Eg, where E is the total
                                   % energy in the system, and randLoss
                                   % is energy lost due tob factors such
                                   % air resistance
        
        if max(Ek) <= 0  % Check if maximum height as been reached, i.e. Ek = 0
            
            E = 0.85 .* abs(Eg);   % No Ek therefore for E = Eg, added 15% energy loss
            
        end
        
        if max(Eg) <= 0  % Check if y = 0, i.e. Eg = m*g*0 = 0;
            
            E = 0.6 .* abs(Ek); % No Eg there for E = Ek, added 40% energy loss due to impact
            Vy = sqrt(2 .* E .* m^-1);  % Calculate new intial velocity from remaining Ek
            
        end
        
        if min(y(c)) < 0  % Check in any projectiles are below ground level
            
            Ylt = y < 0;  % Determine which projectiles are below ground level
            y(Ylt) = 0;  % Set projectiles below ground level to 0
            
        end
        
        if E <= 0.005 .* initialEnergy  % Check if remaining energy in the system is approximately 0          
            break   % END SIMULATION if there is neglegible energy left
        end
        
    end

end

figure(1);
xlim([min(min(x)) max(max(x))]);
ylim([min(min(y)) max(max(y))]);
hold on;

for i = 1 : nTraj
    
    modifiedcomet(x(:,i),y(:,i),0.2,i); % Plot trajectories using the modified comet.m function
    pause(0.2);
    
end

hold off;

% Code can be improved by modelling real energy losses instead of using
% random loss. More than 5 trajectories can be modelled using plot(x,y), or
% by improving the function modifiedcomet.m to allow for more trajectories