function [ output_args ] = PlotVars(c,Limits)
global Vx Vy L W x y Fx Fy C
global Phi nAtoms time Mass0 Mass1 Pty0in Pty1in
global LJEpsilon Phi0 PhiTot KETot MinX MaxX MinY MaxY
global T T0 T1 ScaleV MarkerSize doPlotImage PlotCount
global PlotFig map im PlotSize

if isempty(Limits)
    Limits = Limits;
end


if doPlotImage == 1
    if PlotCount == 0
        close all;
        set(0,'DefaultFigureWindowStyle','normal')
        PlotFig = figure('Position', PlotSize);
        PlotCount = 1;
    else
        PlotCount = PlotCount+1;
    end
end


V2 = Vx.*Vx + Vy.*Vy;
MaxV = max(sqrt(V2));
dx = MaxX - MinX;
dy =  MaxY - MinY;

if c==1
    ScaleV = sqrt(dx*dx + dy*dy)/MaxV*0.15;
end

subplot(3,2,1),plot(x(Pty0in),y(Pty0in),'bo','markers',...
    MarkerSize,'MarkerFaceColor','b');
hold on
subplot(3,2,1),plot(x(Pty1in),y(Pty1in),'go','markers',...
    MarkerSize,'MarkerFaceColor','g');
subplot(3,2,1),quiver(x,y,Fx,Fy,0,'m','linewidth',2);
hold off
axis(Limits);


subplot(3,2,2),plot(time,KETot,'linewidth',2);
hold on
subplot(3,2,2),plot(time,PhiTot+LJEpsilon*nAtoms,'g','linewidth',2);
subplot(3,2,2),plot(time,PhiTot+KETot+LJEpsilon*nAtoms,'k','linewidth',2);
%     axis([0 TStop 0e-22 3e-23]);


subplot(3,2,3),plot(x(Pty0in),y(Pty0in),'bo','markers',...
    MarkerSize,'MarkerFaceColor','b');
hold on
subplot(3,2,3),plot(x(Pty1in),y(Pty1in),'go','markers',...
    MarkerSize,'MarkerFaceColor','g');
subplot(3,2,3),quiver(x,y,Vx*ScaleV,Vy*ScaleV,0,'r','linewidth',2);
hold off
axis(Limits);

AE(Pty0in) = 1/2*Mass0*V2(Pty0in) + Phi(Pty0in)-Phi0;
AE(Pty1in) = 1/2*Mass1*V2(Pty1in) + Phi(Pty1in)-Phi0;
% if AddParticle, AE(end) = 0;end
% AE = Phi;
subplot(3,2,4),scatter3(x,y,AE,ones(1,nAtoms)*300,AE,'fill');
axis(Limits);
view(2);


subplot(3,2,5),plot(x(Pty0in),y(Pty0in),'bo','markers',...
    MarkerSize,'MarkerFaceColor','b');
hold on
subplot(3,2,5),plot(x(Pty1in),y(Pty1in),'go','markers',...
    MarkerSize,'MarkerFaceColor','g');
hold off


subplot(3,2,2),plot(time,KETot,'linewidth',2);
hold on
subplot(3,2,2),plot(time,PhiTot+LJEpsilon*nAtoms,'g','linewidth',2);
subplot(3,2,2),plot(time,PhiTot+KETot+LJEpsilon*nAtoms,'k','linewidth',2);
%     axis([0 TStop 0e-22 3e-23]);


subplot(3,2,6),plot(time,T,'k','linewidth',2);
hold on
subplot(3,2,6),plot(time,T0,'b','linewidth',2);
subplot(3,2,6),plot(time,T1,'g','linewidth',2);
%

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

