function [ U ] = Pot_para(x, paras)

a = paras(1);
x0 = paras(2);
U = paras(1) * (x-x0)^2;

end
