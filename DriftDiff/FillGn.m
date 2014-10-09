function FillGn(nx,dx,SFG)
global Gn n mu Dn V Em DnM muM A B

    muM(1:nx-1) = (mu(1:nx-1) + mu(2:nx))/2;
    DnM(1:nx-1) = (Dn(1:nx-1) + Dn(2:nx))/2;

    Em(1:nx-1) = -(V(2:nx) - V(1:nx-1))/dx;

    Gn = sparse(nx,nx);

if ~SFG
    % J(i) = q n mu E + q Dn dn/dx  % i + 1/2 -- midpoint

    % J(i) = q nM muA EM + q DnM (n(i+1) - n(i))/dx
    % nM = (n(i) + n(i+1))/2

    % J(i) = q ( (n(i) + n(i+1))/2 muA Em + DnM (n(i+1) - n(i))/dx)
    % J(i) = q ( muA Em/2 n(i) + muA Em/2 n(i+1) + DnM/dx n(i+1) - DnM/dx n(i))
    % J(i) = q ((DnM/dx + muA Em/2) n(i+1) - (DnM/dx - muA Em/2) n(i))
    % J(i) = q (A n(i+1) - B n(i))
    % dJn/dx = 0
    % J(i) - J(i-1) = 0
    %  A(i) n(i+1) - B(i) n(i) -  A(i-1) n(i) + B(i-1) n(i-1) = 0

    A = DnM/dx + muM.*Em/2;
    B = DnM/dx - muM.*Em/2;

    for i = 1:nx
        if i == 1
            Gn(i,i) = 1;
        elseif i == nx
            Gn(i,nx-1) = -B(i-1);
            Gn(i,nx) = A(i-1);
        else
            Gn(i,i+1) = A(i);
            Gn(i,i) = -B(i) - A(i-1);
            Gn(i,i-1) = B(i-1);
        end
    end
else

    % J(i) = q n mu E + q Dn dn/dx  % i + 1/2 -- midpoint
    v = muM.*Em;
    vi = v == 0;
    v(vi) = 1e-6;
    % n(i) = A exp(v(i)x/D(i)) + B

    expF = exp(-v.*dx./DnM);
    SCF = 1./(expF - 1);

    % J(i) = v(i) (n(i) expF(i) - n(i+1))*SCF(i)
    % J(i) - J(i-1) = 0
    % v(i) (n(i) expF(i) - n(i+1))*SCF(i) - v(i-1) (n(i-1) expF(i-1) - n(i))*SCF(i) = 0
    % v(i)n(i) expF(i)*SCF(i) - v(i) n(i+1)*SCF(i) - v(i-1)n(i-1)expF(i-1)*SCF + v(i-1)*n(i)*SCF(i) = 0

    % - v(i)*SCF n(i+1) + (v(i-1)+v(i)expF(i))*SCF n(i) - v(i-1)expF(i-1)*SCF(i)n(i-1)= 0


    for i = 1:nx
        if i == 1
            Gn(i,i) = 1;
        elseif i == nx
            Gn(i,nx-1) = v(i-1)*expF(i-1)*SCF(i-1);
            Gn(i,nx) = - v(i-1)*SCF(i-1);
        else
            Gn(i,i+1) = - v(i)*SCF(i);
            Gn(i,i) = (v(i-1)+v(i)*expF(i))*SCF(i);
            Gn(i,i-1) = -v(i-1)*expF(i-1)*SCF(i);
        end
    end

end
end
