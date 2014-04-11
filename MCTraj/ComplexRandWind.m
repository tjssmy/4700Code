function [ Wind ] = ComplexRandWind(nTraj,paras)

% pdf
% wind(x) = 0             x < 0 

% wind(x) = a             x < a

% wind(x) = a - (x-a)    
% wind(x) = 2a - x        x < 2a

% wind(x) = 0             x < inf

% cdf
% c = a*x                                0 < c < a^2
% c = a^2 - 2a^2 - 4a^2/2 + 2ax + x^2/2  
% c = -3a^2 + 2ax + x^2/2                a^2 < c < 3a^2
% c = a^2                                c > 3a^2

a = paras(1);
M = 3*a^2;
r = rand(1,nTraj)*M;


ri = r < a^2;
Wind(ri) = r(ri)/a;

hist(Wind);

rin = ~ri;
%  0 = -(c+3a^2) + 2ax + x^2/2                a^2 < c < 3a^2
Wind(rin) = -2*a + sqrt(4*a^2 + 4/2*(r(rin)+3*a^2));

hist(Wind);

end

