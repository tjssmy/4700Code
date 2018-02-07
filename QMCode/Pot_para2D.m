function [ U ] = Pot_para2D(x,y, paras)

a = paras(1);
x0 = paras(2);
y0 = paras(3);
U = a * ((x-x0)^2 + (y-y0)^2);

end
