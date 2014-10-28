function [ U ] = Pot_step(x, paras)

dx = paras(1);
x0 = paras(2);
a = paras(3);
b = paras(4);

if (x < x0-dx)
    U = a;
elseif (x > x0 + dx)
    U = a;
else
    U = b;
end

end
