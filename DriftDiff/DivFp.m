function DivFp(nx,dx,BC)
global pp Em DpM MupM divFp
    
    
    % J(i) = q p Mup E - q Dp dp/dx  % i + 1/2 -- midpoint
    
    % J(i) = q pM MupA EM - q DpM (p(i+1) - p(i))/dx
    % pM = (p(i) + p(i+1))/2
    
    % J(i) = q ( (p(i) + p(i+1))/2 MupA Em - DpM (p(i+1) - p(i))/dx)
    % J(i) = q ( MupA Em/2 p(i) + MupA Em/2 p(i+1) - DpM/dx p(i+1) + DpM/dx p(i))
    % J(i) = q ((-DpM/dx + MupA Em/2) p(i+1) - (-DpM/dx - MupA Em/2) p(i))
    % J(i) = q (A p(i+1) - B p(i))
    % dJp/dx = (J(i) - J(i-1))/dx 
    
    %  dJp/dx = q/dx((A(i) p(i+1) - B(i) p(i) -  A(i-1) p(i) + B(i-1) p(i-1)) 
    %  dFp/dx = ((A(i) p(i+1) - B(i) p(i) -  A(i-1) p(i) + B(i-1) p(i-1))/dx 
    
    A = -DpM/dx + MupM.*Em/2;
    B = -DpM/dx - MupM.*Em/2;
    
    for i = 1:nx
        if i == 1
            if BC == 0
                divFp(i) = (A(i)*pp(i+1) - B(i)*pp(i) - 0)/dx;
            else
                divFp(i) = 0;
            end
        elseif i == nx
            if BC == 0
            divFp(i) = (0 - A(i-1)*pp(i) + B(i-1)*pp(i-1))/dx;
            else
                divFp(i) = 0;
            end
        else
            divFp(i) = (A(i)*pp(i+1) - B(i)*pp(i) - ...
                A(i-1)*pp(i) + B(i-1)*pp(i-1))/dx;
        end
    end
    
    divFp = -divFp;
end



