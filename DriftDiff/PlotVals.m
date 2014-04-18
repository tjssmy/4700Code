function PlotVals(nx,dx,h,axis)
global V Em n p x xm C DnM MunM DpM MupM divFn divFp Rho niSi t
global TwoCarriers


gradn = (n(2:nx)-n(1:nx-1))/dx;
nM = (n(2:nx)+n(1:nx-1))/2;

JnDiff = C.q_0*DnM.*gradn;
JnDrift = C.q_0*MunM.*nM.*Em;

gradp = (p(2:nx)-p(1:nx-1))/dx;
pM = (p(2:nx)+p(1:nx-1))/2;

JpDiff = C.q_0*DpM.*gradp;
JpDrift = C.q_0*MupM.*pM.*Em;


subplot(4,3,1),plot(x,V);
title('V');
%     axis([0 l -0.05 0.01])
hold(h)
subplot(4,3,2),plot(xm,Em);
title('E');
%     axis([0 l -15e4 1e4])
hold(h)

subplot(4,3,3),plot(x,Rho);
title('Rho');
%     axis([0 l -15e4 1e4])
hold(h)

subplot(4,3,4),plot(x,n);
title('n');
hold(h)

subplot(4,3,5),plot(x,p);
title('p');
hold(h)

subplot(4,3,6),plot(x,n.*p - niSi^2);
title('np-ni^2');
hold(h)

subplot(4,3,7),plot(x,divFn);
title('divFn');
hold(h)

subplot(4,3,8),plot(x,divFp);
title('divFp');
hold(h)


%     axis([0 l 0 5e20])
subplot(4,3,9),plot(xm,-JnDiff,'r');
hold on
subplot(4,3,9),plot(xm,JnDrift,'b');
title('JnDiff/-JnDrift');
%     axis([0 l -5e5 5e5])
hold(h);    

%     axis([0 l 0 5e20])
subplot(4,3,10),plot(xm,JpDiff,'r');
hold on
subplot(4,3,10),plot(xm,JpDrift,'b');
title('JpDiff/-JpDrift');
%     axis([0 l -5e5 5e5])
hold(h);    

Jerr = (JnDrift+JnDiff)/max(abs([JnDrift, JnDiff])*100);
subplot(4,3,11),plot(xm,Jerr);
hold on
if TwoCarriers ==1 
    Jerr = (JpDrift+JpDiff)/max(abs([JpDrift, JpDiff])*100);
    subplot(4,3,11),plot(xm,Jerr,'r');
end
title('JDiff+JDrift (normalized %)');
%     axis([0 l -5e5 5e5])
hold(h);    

subplot(4,3,12),plot(t,max(n),'+b');
hold on
subplot(4,3,12),plot(t,max(p),'or');
title('Max n and p');
hold on

end

