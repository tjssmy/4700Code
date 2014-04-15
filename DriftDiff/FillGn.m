function FillGn(nx,dx)
global Gn n mu Dn V


% J(i) = - q n mu E + q Dn dn/dx  % i + 1/2 -- midpoint

% J(i) = - q nM muA EM + q DnM (n(i+1) - n(i))/dx
% nM = (n(i) + n(i+1))/2

% J(i) = q ( - (n(i) + n(i+1))/2 muA Em + DnM (n(i+1) - n(i))/dx)
% J(i) = q ( - muA Em/2 n(i) - muA Em/2 n(i+1) + DnM/dx n(i+1) - DnM/dx n(i))
% J(i) = q ((DnM/dx - muA Em/2) n(i+1) - (muA Em/2 + DnM/dx) n(i))
% J(i) = q (A n(i+1) - B n(i))
% dJn/dx = 0
% J(i) - J(i-1) = 0
%  A(i) n(i+1) - B(i) n(i) -  A(i-1) n(i) + B(i-1) n(i-1) = 0

muA(1:nx-1) = (mu(1:nx-1) + mu(2:nx))/2
DnM(1:nx-1) = (Dn(1:nx-1) + Dn(2:nx))/2
EM(1:nx-1) = -(V(1:nx-1) + V(2:nx))/dx;

A = DnM/dx - muA.*EM/2;
B = muA.*EM/2 + DnM/dx;

for i = 1:nx
    if i == 1              % J(1) = 0;
        Gn(i,i) = -B(i);
        Gn(i,i+1) = A(i);
    elseif i == nx         %(J(nx) = 0
        Gn(i,nx-1) = -B(i-1);
        Gn(i,nx) = A(i-1);
    else
        Gn(i,i+1) = A(i-1);
        Gn(i,i) = -B(i) - A(i-1);
        Gn(i,i-1) = B(i-1);
    end
end
end



