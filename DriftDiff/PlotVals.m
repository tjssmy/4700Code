function PlotVals(nx,dx,h,l,TStop,Yaxis)
global V Em n p x xm C DnM MunM DpM MupM divFn divFp Rho niSi t
global TwoCarriers PlotSS PlotFig PlotCount im map doPlotImage

if doPlotImage == 1
    if PlotCount == 0
        close all;
        set(0,'DefaultFigureWindowStyle','normal')
        PlotFig = figure('Position', [100, 100, 1049, 895]);
        PlotCount = 1;
    else
        PlotCount = PlotCount+1;
    end
end
    
if isempty(Yaxis)
    scale = 0;
else
    scale = 1;
end
    
gradn = (n(2:nx)-n(1:nx-1))/dx;
nM = (n(2:nx)+n(1:nx-1))/2;

JnDiff = C.q_0*DnM.*gradn;
JnDrift = C.q_0*MunM.*nM.*Em;

gradp = (p(2:nx)-p(1:nx-1))/dx;
pM = (p(2:nx)+p(1:nx-1))/2;

JpDiff = -C.q_0*DpM.*gradp;
JpDrift = C.q_0*MupM.*pM.*Em;


subplot(4,3,1),plot(x,V);
title('V');
if scale, axis([ 0 l Yaxis{1}]); end; 

hold(h)
subplot(4,3,2),plot(xm,Em);
if scale, axis([ 0 l Yaxis{2}]); end
title('E');
hold(h)

subplot(4,3,3),plot(x,Rho);
title('Rho');
if scale, axis([ 0 l Yaxis{3}]); end
hold(h)

subplot(4,3,4),plot(x,n);
if scale, axis([ 0 l Yaxis{4}]); end
title('n');
hold(h)

subplot(4,3,5),plot(x,p);
if scale, axis([ 0 l Yaxis{5}]); end
title('p');
hold(h)

subplot(4,3,6),plot(x,n.*p - niSi^2);
if scale, axis([ 0 l Yaxis{6}]); end
title('np-ni^2');
hold(h)

subplot(4,3,7),plot(x,divFn);
if scale, axis([ 0 l Yaxis{7}]); end
title('divFn');
hold(h)

subplot(4,3,8),plot(x,divFp);
if scale, axis([ 0 l Yaxis{8}]); end
title('divFp');
hold(h)

subplot(4,3,9),plot(xm,-JnDiff,'r');
hold on
subplot(4,3,9),plot(xm,JnDrift,'b');
if scale, axis([ 0 l Yaxis{9}]); end
title('JnDiff/-JnDrift');
hold(h);    

subplot(4,3,10),plot(xm,-JpDiff,'r');
hold on
subplot(4,3,10),plot(xm,JpDrift,'b');
if scale, axis([ 0 l Yaxis{10}]); end
title('JpDiff/-JpDrift');
hold(h);    

if PlotSS
Jerr = (JnDrift+JnDiff)/max(abs([JnDrift, JnDiff])*100);
subplot(4,3,11),plot(xm,Jerr);
hold on
if TwoCarriers ==1 
    Jerr = (JpDrift+JpDiff)/max(abs([JpDrift, JpDiff])*100);
    subplot(4,3,11),plot(xm,Jerr,'r');
end
title('JDiff+JDrift (normalized %)');
if scale, axis([ 0 l Yaxis{11}]); end
hold(h);    
else
    if TwoCarriers ==1
        Jtot = JpDrift+JpDiff+JnDrift+JnDiff;
    else
        Jtot = JnDrift+JnDiff;
    end
    subplot(4,3,11),plot(xm,Jtot,'r');
    title('Jtot');
    if scale, axis([ 0 l Yaxis{11}]); end
    hold(h);
    
end
subplot(4,3,12),plot(t,max(n),'+b');
hold on
subplot(4,3,12),plot(t,max(p),'or');
if scale, axis([ 0 TStop Yaxis{12}]); end
title('Max n and p');
hold on

if doPlotImage   
    if PlotCount == 1
        f = getframe(PlotFig);
        [im,map] = rgb2ind(f.cdata,256,'nodither');
        im(1,1,1,2) = 0;
    else
        f = getframe(PlotFig);
        im(:,:,1,PlotCount) = rgb2ind(f.cdata,map,'nodither');
    end
end

end

