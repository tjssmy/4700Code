% clear all
set(0,'DefaultFigureWindowStyle','docked')
global C

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light

% - hb^2/2m dP^2/dx^2 + V(x)P = E P

nx = 200;
l = .40e-9;
x =linspace(0, l, nx);
dx = x(2) - x(1);

pot = @Pot_const;
paras = [0 * C.q_0];

% pot = @Pot_well;
% paras = [l/4,l/2,75*C.q_0,0*C.q_0];

% pot = @Pot_Dwell;
% paras = [l/16,l/4,l/2,75*C.q_0,0*C.q_0];

% pot = @Pot_para;
% paras = [1000/(l^2)*C.q_0,l/2];

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

[V,D] = eig(G);

[Ds,Pr] = sort(diag(D));
D = D(Pr,Pr);
V = V(:,Pr);

subplot(2,1,1),plot(x*1e9,Poten/C.q_0);
xlabel('x (nm)');
ylabel('Potential (eV)');

subplot(2, 1, 2), plot(x * 1e9, V(:, 1:5))
xlabel('x (nm)');
ylabel('Eigenvectors (n/a)');

E(1) = C.hb^2 * pi * pi / (2 * C.m_0 * l^2) / C.q_0;
E(2) = 4 * C.hb^2 * pi * pi / (2 * C.m_0 * l^2) / C.q_0;
E(3) = 9 * C.hb^2 * pi * pi / (2 * C.m_0 * l^2) / C.q_0;

display(E)
display((diag(D(1:8, 1:8)) / C.q_0)')
