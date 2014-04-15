function [ Wind ] = ComplexRandWind(nTraj,paras)

% pdf
% wind(x) = 0               x < 0 

% wind(x) = 1.0             x < a

% wind(x) = 1.0 - 1/a*(x-a) 
% wind(x) = 2.0 - 1/a*x     x < 2a   

% wind(x) = 0               x > 2a

% cdf
% c = x                                           0 < c < a
% c = a + 2(x) - 1/(2a)*(x)^2 - 2.0a + 1/(2a)*a^2  
% c =  -1/2a + 2(x) - 1/(2a)*(x)^2                a < c < 1.5a
% c = 1.5a                                        c > 1.5a

a = paras(1);
M = 1.5*a;
r = rand(1,nTraj)*M;

ri = r < a;
Wind(ri) = r(ri);

rin = ~ri;
% -(c+1/2a) + 2(x) - 1/(2a)*(x)^2 = 0             
A = -1/(2*a);
B = 2;
C = -(r(rin)+a/2);
Wind(rin) = (-B + sqrt(B*B - 4*A*C))/(2*A);


% hist(Wind,50);
% xlabel('Wind');
% ylabel('number');
Wind = -Wind;

end

