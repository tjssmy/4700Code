function [ U ] = Pot_Dwell(x, paras)

dx = paras(1);
dw = paras(2);
x0 = paras(3);
a = paras(4);
b = paras(5);

if (x < x0-dw-dx/2)
    U = a;
elseif (x < x0-dx/2)
    U = b;
elseif (x < x0+dx/2)
    U = a;
elseif (x < x0+dx/2+dw)
    U = b;
else 
    U = a;
end

end
