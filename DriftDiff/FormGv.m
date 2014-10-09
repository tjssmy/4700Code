function  FormGv(nx,lBC,rBC)
global Gv Bv

% Poisson equation d^2V/dx^2 = -1/EpiSi rho(x)
% Gv V = -dx^2/EpiSi rho(x)
% E = - dV/dx

Gv = sparse(nx,nx);
Bv = zeros(1,nx);

for i = 1:nx

    if i == 1
        if lBC == 'fl'
            Gv(i,i) = 1;
            Gv(i,i+1) = -1; % floating contact
        else
            Bv(i) = lBC;
            Gv(i,i) = 1;
        end
    elseif i == nx
        Bv(i) = rBC;
        Gv(i,i) = 1;  % fixed contact RVbc
    else
        Gv(i,i-1) = 1;
        Gv(i,i) = -2;
        Gv(i,i+1) = 1;
    end
end

end
