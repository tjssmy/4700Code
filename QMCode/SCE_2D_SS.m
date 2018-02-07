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

nx = 50;
l = .40e-9;
x =linspace(0, l, nx);

d = x(2) - x(1);

w = .50e-9;
ny = round(w/d)+1;
w = (ny-1)*d;
y =linspace(0, w, ny);


% pot = @Pot_const2D;
% paras = [0 * C.q_0];

pot = @Pot_well2D;
paras = [l/4,l/2,w/4,w/2,175*C.q_0,0*C.q_0];

% pot = @Pot_Dwell;
% paras = [l/16,l/4,l/2,75*C.q_0,0*C.q_0];

% pot = @Pot_para;
% paras = [1000/(l^2)*C.q_0,l/2];

d2 = d^2;
B = (C.hb^2) / (2 * C.m_0);
G = sparse(nx*ny,nx*ny);

for i = 1:nx
    for j = 1:ny
        n = i + (j - 1) * nx;
        
        nxm = i - 1 + (j-1)*nx;
        nxp = i + 1 + (j-1)*nx;
        nym = i + (j-1-1)*nx;
        nyp = i + (j-1+1)*nx;
        
        Poten(i,j) = pot(x(i), y(j), paras);
               
        if i == 1
            G(n, n) = B / d2;
        elseif i == nx
            G(n, n) = B / d2;
        elseif j == 1
            G(n, n) = B / d2;
        elseif j == ny
            G(n, n) = B / d2;
        else
            G(n, n) = B * 4 / d2 + Poten(i,j);
            G(n,nxm) = -B / d2;
            G(n,nxp) = -B / d2;
            G(n,nym) = -B / d2;
            G(n,nyp) = -B / d2;
        end
        
    end
end

[V,D] = eigs(G,8,0);

% [Ds,Pr] = sort(diag(D));
% D = D(Pr,Pr);
% V = V(:,Pr);

Mmap = zeros(nx,ny,8);
for m = 1:8
%     for i = 1:nx
%         for j = 1:ny
%             n = i + (j - 1) * nx;
%             Mmap(i, j, m) = V(n,m);
%         end
%     end
    Mmap(:,:,m) = reshape(V(:,m),[nx,ny]);
end

subplot(3,3,1),surf(x*1e9,y*1e9,Poten.'/C.q_0,'linestyle','none');
xlabel('x (nm)');
xlabel('y (nm)');
zlabel('Potential (eV)');


for m = 1:8
    subplot(3, 3, m+1), surf(x * 1e9,y*1e9,real(Mmap(:,:,m)).','linestyle','none')
    xlabel('x (nm)');
    xlabel('y (nm)');
    zlabel('\psi');
    title(num2str(D(m,m)))
end

