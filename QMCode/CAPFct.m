function [ F ] = CAPFct(x, paras)

del = paras(1);
A = paras(2);

F(1:length(x)) = 0;

ii = x < del;
F(ii) = -1i * A * (del - x(ii)).^2;
ii = x > x(end) - del;
F(ii) = -1i * A * ((x(end) - del) - x(ii)).^2;

end
