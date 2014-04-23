function [ Phi dPhidr ] = LJPot(r,Epsilon,Sigma)

Phi = 4*Epsilon*(Sigma^12*r^-12 - Sigma^6*r^-6);

dPhidr = 4*Epsilon*(Sigma^12*(-12)*r^-13 - 14Sigma^6*(-6)*r^-7);

end

