function PlotVals(nx,dx,h,l,TStop,Yaxis,pl)
global V Em n p x xm C DnM MunM DpM MupM divFn divFp Rho niSi t
global TwoCarriers PlotSS PlotFig PlotCount im map doPlotImage Itot
global tv V1 JnDiff JpDiff JnDrift JpDrift Jtot


gradn = (n(2:nx)-n(1:nx-1))/dx;
nM = (n(2:nx)+n(1:nx-1))/2;

JnDiff = C.q_0*DnM.*gradn;
JnDrift = C.q_0*MunM.*nM.*Em;

gradp = (p(2:nx)-p(1:nx-1))/dx;
pM = (p(2:nx)+p(1:nx-1))/2;

JpDiff = -C.q_0*DpM.*gradp;
JpDrift = C.q_0*MupM.*pM.*Em;


if TwoCarriers ==1
    Jtot = JpDrift+JpDiff+JnDrift+JnDiff;
else
    Jtot = JnDrift+JnDiff;
end

Itot = [Itot [Jtot(1); Jtot(end)]];
V1 = [V1 V(1)];
tv = [tv t];

if ~pl
    return
end

if isempty(Yaxis)
    scale = 0;
else
    scale = 1;
end


subplot(3,3,1),plot(x,V);
title('V');
if scale, axis([ 0 l Yaxis{1}]); end;

hold(h)
subplot(3,3,2),plot(xm,Em);
if scale, axis([ 0 l Yaxis{2}]); end
title('E');
hold(h)

subplot(3,3,3),plot(x,Rho);
title('Rho');
if scale, axis([ 0 l Yaxis{3}]); end
hold(h)

subplot(3,3,4),plot(x,n);
if scale, axis([ 0 l Yaxis{4}]); end
title('n');
hold(h)

subplot(3,3,5),plot(x,p);
if scale, axis([ 0 l Yaxis{5}]); end
title('p');
hold(h)

if TwoCarriers ==1
    Jtot = JpDrift+JpDiff+JnDrift+JnDiff;
else
    Jtot = JnDrift+JnDiff;
end
subplot(3,3,6),plot(xm,Jtot,'r');
title('Jtot');
if scale, axis([ 0 l Yaxis{11}]); end
hold(h);

Itot = [Itot [Jtot(1); Jtot(end)]];
V1 = [V1 V(1)];
tv = [tv t];

subplot(3,3,7)
% plot(tv,log10(abs(Itot(1,:)))); hold on
% plot(tv,log10(abs(Itot(2,:)))); hold off
plot(tv,-Itot(1,:)); hold on
plot(tv,-Itot(2,:)); hold off

title('Itot')

subplot(3,3,8)
plot(tv,V1); 
title('V1')


end


