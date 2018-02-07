function [ U ] = Pot_well2D(x,y, paras)

dx = paras(1);
x0 = paras(2);
dy = paras(3);
y0 = paras(4);

a = paras(5);
b = paras(6);

if x > x0-dx && x < x0 + dx && y > y0 - dy && y < y0 + dy
    U = b;
else
    U = a;
end

end
