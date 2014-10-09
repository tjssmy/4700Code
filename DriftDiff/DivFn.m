function DivFn(nx, dx, BC)
global np Em DnM MunM divFn

    % J(i) = q n Mun E + q Dn dn / dx  % i + 1/2 -- midpoint

    % J(i) = q nM MunA EM + q DnM (n(i + 1) - n(i)) / dx
    % nM = (n(i) + n(i + 1))/2

    % J(i) = q ( (n(i) + n(i + 1))/2 MunA Em + DnM (n(i + 1) - n(i)) / dx)
    % J(i) = q ( MunA Em / 2 n(i) + MunA Em / 2 n(i + 1) + DnM / dx n(i + 1) - DnM / dx n(i))
    % J(i) = q ((DnM / dx + MunA Em / 2) n(i + 1) - (DnM / dx - MunA Em / 2) n(i))
    % J(i) = q (A n(i + 1) - B n(i))
    % dJn / dx = (J(i) - J(i - 1)) / dx

    %  dJn / dx = q / dx((A(i) n(i + 1) - B(i) n(i) -  A(i - 1) n(i) + B(i - 1) n(i - 1))
    %  dFn / dx = ((A(i) n(i + 1) - B(i) n(i) -  A(i - 1) n(i) + B(i - 1) n(i - 1)) / dx

    A = DnM / dx + MunM .* Em / 2;
    B = DnM / dx - MunM .* Em / 2;

    for i = 1:nx
        if i == 1
            if BC == 0
                divFn(i) = (A(i) * np(i + 1) - B(i) * np(i) - 0) / dx;
            else
                divFn(i) = 0;
            end
        elseif i == nx
            if BC == 0
                divFn(i) = (0 - A(i - 1) * np(i) + B(i - 1) * np(i - 1)) / dx;
            else
                divFn(i) = 0;
            end
        else
            divFn(i) = (A(i) * np(i + 1) - B(i) * np(i) - ...
                A(i - 1) * np(i) + B(i - 1) * np(i - 1)) / dx;
        end
    end
end
