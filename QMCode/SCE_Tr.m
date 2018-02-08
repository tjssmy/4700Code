clearvars
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


nx = 501;
len = 10e-9;
x =linspace(0, len, nx);
dx = x(2) - x(1);

pot = @Pot_const;
paras = [0*C.q_0];

% pot = @Pot_well;
% paras = [len/64,len/2,0*C.q_0,1*C.q_0];

% pot = @Pot_Dwell;
% paras = [len / 16, len / 64, len / 2, 0 * C.q_0, 10 * C.q_0];

pot = @Pot_para;
paras = [10/(len^2)*C.q_0,len/2];




dx2 = dx^2;
B = (C.hb^2) / (2 * C.m_0);
G = sparse(nx,nx);

% Cap = CAPFct(x, [len * 0.15, 1]);
Cap(1:nx) = 0;

for i = 1:nx
%     Poten(i) = pot(x(i),paras);
    Poten(i) = pot(x(i), paras) + Cap(i);

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
% subplot(2,1,1), plot(x,Poten/C.q_0);
% subplot(2,1,2), plot(x,V(:,1:5))
%
% figure

% Gi = inv(G);

subplot(3, 1, 1), plot(x, real(Poten / C.q_0), 'k');
hold
subplot(3, 1, 1), plot(x, imag(Poten / C.q_0), 'r');

nSteps = 300;
Tmax = 5e-15;
Tmax0 = 1e-13;


K0 = C.m_0 * (len / Tmax0) / C.hb * 20; % P = hb*k k = mv/hb;
% K0 = 0;
P0(1, 1:nx) = (exp(-(x - len / 4).^2 / ((len / 50)^2))) .* exp(1i * K0 * x);
Pp = P0.';
maxP0 = max(abs(P0));
subplot(3, 1, 2), plot(x, real(Pp), 'r');
hold on
subplot(3, 1, 2), plot(x, imag(Pp), 'b');
subplot(3, 1, 2), plot(x, abs(Pp), 'k');
axis([0 len -maxP0 maxP0])
hold off

f = getframe;
[im,map] = rgb2ind(f.cdata, 256, 'nodither');
im(1, 1, 1, 20) = 0;

% set(gca,'nextplot','replacechildren','visible','off')

Norm=zeros(1, nSteps);
Norm(1) = sum(abs(Pp).^2);
subplot(3, 1, 3), plot(Norm, 'b');
axis([0 nSteps 0 Norm(1)*1.1])
dt = Tmax / nSteps;

H = speye(nx) - 1i * dt / C.hb * G;
Hp = speye(nx) + 1i * dt / C.hb * G;

[l,u] = lu(Hp);

for j = 1:nSteps
    %      P = u\(l\(Pp'));
    P = u\(l\(H*Pp));

    subplot(3, 1, 2), plot(x, real(P), 'r');
    hold on
    subplot(3, 1, 2), plot(x, imag(P), 'b');
    subplot(3, 1, 2), plot(x, abs(P), 'k');
    axis([0 len -maxP0 maxP0])
    hold off

    f = getframe;
    im(:, :, 1, j) = rgb2ind(f.cdata, map, 'nodither');

    Norm(j) = sum(abs(Pp).^2);
    subplot(3, 1, 3), plot(Norm, 'b');
    axis([0 nSteps 0 Norm(1)*1.1])

    %     pause(0.004)

    Pp = P;
end

imwrite(im, map, 'imagefile.gif', 'DelayTime', 0, 'LoopCount', inf);

% subplot(3, 1, 2), plot(x,V(:,1:5))
