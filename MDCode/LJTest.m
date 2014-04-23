clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle','docked')

AtomSpacing = 0.5430710e-9;
LJSigma = AtomSpacing/2^(1/6);
LJEpsilon = 1e-21;

for n = 150:1000
   
    r = n*LJSigma/200;
    [Phi dPhidr] = LJPot(r,LJEpsilon,LJSigma);
    
    PhiT(n) = Phi;
    dPhidrT(n) = dPhidr;
    rT(n) = r;
end

subplot(2,1,1),plot(rT,PhiT,'x');
subplot(2,1,2),plot(rT,-dPhidrT,'x');