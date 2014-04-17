function [divJn] = DivJn(nx,dx)
global Gn np mu Dn V Em DnM muM A B
    
    muM(1:nx-1) = (mu(1:nx-1) + mu(2:nx))/2;
    DnM(1:nx-1) = (Dn(1:nx-1) + Dn(2:nx))/2;
    
    Em(1:nx-1) = -(V(2:nx) - V(1:nx-1))/dx;


    % J(i) = q n mu E + q Dn dn/dx  % i + 1/2 -- midpoint
    
    % J(i) = q nM muA EM + q DnM (n(i+1) - n(i))/dx
    % nM = (n(i) + n(i+1))/2
    
    % J(i) = q ( (n(i) + n(i+1))/2 muA Em + DnM (n(i+1) - n(i))/dx)
    % J(i) = q ( muA Em/2 n(i) + muA Em/2 n(i+1) + DnM/dx n(i+1) - DnM/dx n(i))
    % J(i) = q ((DnM/dx + muA Em/2) n(i+1) - (DnM/dx - muA Em/2) n(i))
    % J(i) = q (A n(i+1) - B n(i))
    % dJn/dx = (J(i) - J(i-1))/dx 
    %  dJn/dx = q*dx((A(i) n(i+1) - B(i) n(i) -  A(i-1) n(i) + B(i-1) n(i-1)) 
    
    A = DnM/dx + muM.*Em/2;
    B = DnM/dx - muM.*Em/2;
    
    divJn = zeros(1,nx);
    for i = 1:nx
        if i == 1
            divJn(i) = q*dx*(A(i)*np(i+1) - B(i)*np(i) - 0);
        elseif i == nx
            divJn(i) = q*dx*(0 - A(i-1)*np(i) + B(i-1)*np(i-1));
        else
            divJn(i) = q*dx*(A(i)*np(i+1) - B(i)*np(i) - ...
                A(i-1)*np(i) + B(i-1)*np(i-1));
        end
    end
    
end



