% clear all
close all
% set(0,'DefaultFigureWindowStyle','docked')
global C

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light

% -i hb dP/dt = - hb^2/2m dP^2/dx^2 + V(x)P
% -i hb (P - Pp)/dt = G P
% BE:
%   (P - Pp) =  i dt/hb G P
%   (1-i dt/hb G)*P =   Pp
%   H*P =   Pp
%   P = inv(H)*Pp


nx = 2801;
len = 60e-9;
x = linspace(0, len, nx);
dx = x(2) - x(1);

% pot = @Pot_const;
% paras = [0*C.q_0];

pot = @Pot_well;
paras = [len / 64, len / 2, 0 * C.q_0, 1.2 * C.q_0];

% pot = @Pot_Dwell;
% paras = [len/16,len/64,len/2,0*C.q_0,1*C.q_0];

% pot = @Pot_para;
% paras = [10/(len^2)*C.q_0,len/2];

dx2 = dx^2;
B = (C.hb^2) / (2 * C.m_0);

for i = 1:nx
    Poten(i) = pot(x(i), paras);
    if i == 1
        G(i, i) = B / dx2;
    elseif i == nx
        G(i, i) = B / dx2;
    else
        G(i, i) = B * 2 / dx2 + Poten(i);
        G(i, i - 1) = -B / dx2;
        G(i, i + 1) = -B / dx2;
    end

end

% [V,D] = eig(G);
%
% [Ds,Pr]=sort(diag(D));
% D=D(Pr,Pr);
% V=V(:,Pr);
%
% subplot(2,1,1),plot(x,Poten/C.q_0);
% subplot(2,1,2),plot(x,V(:,1:5))
%
% figure

subplot(3, 1, 1), plot(x, Poten / C.q_0, 'k');

nSteps = 200;
Tmax = 1e-13;

K0 = C.m_0 * (len / Tmax) / C.hb * 1; % P = hb*k k = mv/hb;
% K0 = 0;
P0(1, 1:nx) = (exp(-(x - len / 4).^2 / ((len / 32)^2))) .* exp(-1i * K0 * x);
Pp = P0;

Norm=zeros(1, nSteps);
Norm(1) = sum(abs(Pp).^2);

H = eye(nx) - 1i * dt / C.hb * G;
[l, u] = lu(H);

Hp = eye(nx) + 1i * dt / C.hb * G;
[lp, up] = lu(Hp);

for j = 1:nSteps
    %      P = u\(l\(Pp')); % Backward Euler
    P = up\(lp\(H*Pp'));      % Crank-Nickolson
    Pp = P';

    Norm(j) = sum(abs(Pp).^2);
end

imwrite(im, map, 'imagefile.gif', 'DelayTime', 0, 'LoopCount', inf);

% subplot(3,1,2),plot(x,V(:,1:5))
