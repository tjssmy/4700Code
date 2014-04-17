function PlotVals(nx,dx,h,axis)
global V Em n x xm C DnM muM divFn Rho


gradn = (n(2:nx)-n(1:nx-1))/dx;
nM = (n(2:nx)+n(1:nx-1))/2;

JnDiff = C.q_0*DnM.*gradn;
JnDrift = C.q_0*muM.*nM.*Em;


subplot(3,3,1),plot(x,V);
title('V');
%     axis([0 l -0.05 0.01])
hold(h)
subplot(3,3,2),plot(xm,Em);
title('E');
%     axis([0 l -15e4 1e4])
hold(h)

subplot(3,3,3),plot(x,Rho);
title('Rho');
%     axis([0 l -15e4 1e4])
hold(h)

subplot(3,3,4),plot(x,n);
title('n');
hold(h)

subplot(3,3,7),plot(x,divFn);
title('divFn');
hold(h)

%     axis([0 l 0 5e20])
subplot(3,3,8),plot(xm,JnDiff,'r');
hold on
subplot(3,3,8),plot(xm,JnDrift,'b');
title('JnDiff/JnDrift');
%     axis([0 l -5e5 5e5])
hold(h);    

Jerr = (JnDrift+JnDiff)/max(abs([JnDrift, JnDiff])*100);
subplot(3,3,9),plot(xm,Jerr);
title('JnDiff+JnDrift (normalized %)');
%     axis([0 l -5e5 5e5])
hold(h);    

end

