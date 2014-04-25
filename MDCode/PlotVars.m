function [ output_args ] = PlotVars(c)
global Vx Vy L W x y Fx Fy 
global Phi nAtoms time Mass0 Mass1 Pty0in Pty1in
global LJEpsilon Phi0 PhiTot KETot MinX MaxX MinY MaxY

V2 = Vx.*Vx + Vy.*Vy;
MaxV = max(sqrt(V2));
dx = MaxX - MinX;
dy =  MaxY - MinY;
ScaleV = sqrt(dx*dx + dy*dy)/MaxV*0.5;

subplot(2,2,1),plot(x(Pty0in),y(Pty0in),'bo','markers',20,'MarkerFaceColor','b');
hold on
subplot(2,2,1),plot(x(Pty1in),y(Pty1in),'go','markers',20,'MarkerFaceColor','g');
subplot(2,2,1),quiver(x,y,Fx,Fy,0,'m','linewidth',2);
hold off
axis([MinX,MaxX,MinY,MaxY]);


subplot(2,2,3),plot(x(Pty0in),y(Pty0in),'bo','markers',20,'MarkerFaceColor','b');
hold on
subplot(2,2,3),plot(x(Pty1in),y(Pty1in),'go','markers',20,'MarkerFaceColor','g');
subplot(2,2,3),quiver(x,y,Vx*ScaleV,Vy*ScaleV,0,'r','linewidth',2);
hold off
axis([MinX,MaxX,MinY,MaxY]);

subplot(2,2,2),plot(time,KETot,'linewidth',2);
hold on
subplot(2,2,2),plot(time,PhiTot+LJEpsilon*nAtoms,'g','linewidth',2);
subplot(2,2,2),plot(time,PhiTot+KETot+LJEpsilon*nAtoms,'k','linewidth',2);
%     axis([0 TStop 0e-22 3e-23]);

AE = 1/2*Mass0*V2 + Phi-Phi0;
% if AddParticle, AE(end) = 0;end
% AE = Phi;
subplot(2,2,4),scatter3(x,y,AE,ones(1,nAtoms)*300,AE,'fill');
axis([MinX,MaxX,MinY,MaxY]);
view(2);

end

