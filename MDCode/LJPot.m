function [ Phi dPhidr ] = LJPot(r,alpha,rm)

Phi = alpha*((rm/r)^12 - 2*(rm/r)^6);
dPhidr = alpha*(-12*rm*(rm/r)^13 - 2*(-6*rm)*(rm/r)^7);

end

