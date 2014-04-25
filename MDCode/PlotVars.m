function [ output_args ] = PlotVars(c, AddParticle)
global Vx Vy L W x y Fx Fy AtomSpacing LAtoms WAtoms
global Phi nAtoms time Mass0 Mass1 Pty0in Pty1in
global LJEpsilon Phi0 PhiTot KETot

V2 = Vx.*Vx + Vy.*Vy;
MaxV = max(sqrt(V2));
ScaleV = sqrt(L*L + W*W)/MaxV*0.5;

subplot(2,2,1),plot(x(Pty0in),y(Pty0in),'bo','markers',20,'MarkerFaceColor','b');
hold on
subplot(2,2,1),plot(x(Pty1in),y(Pty1in),'go','markers',20,'MarkerFaceColor','g');
subplot(2,2,1),quiver(x,y,Fx,Fy,0,'m','linewidth',2);
hold off

axis([-2*AtomSpacing,(LAtoms+2)*AtomSpacing,...
    -2*AtomSpacing,(WAtoms+2)*AtomSpacing]);


PhiTot(c) = sum(Phi)/2;

KETot(c) = 1/2*Mass0*...
    sum(Vx(Pty0in).*Vx(Pty0in)+Vy(Pty0in).*Vy(Pty0in))...
    + 1/2*Mass1*...
    sum(Vx(Pty1in).*Vx(Pty1in)+Vy(Pty1in).*Vy(Pty1in));

subplot(2,2,3),plot(x(Pty0in),y(Pty0in),'bo','markers',20,'MarkerFaceColor','b');
hold on
subplot(2,2,3),plot(x(Pty1in),y(Pty1in),'go','markers',20,'MarkerFaceColor','g');
subplot(2,2,3),quiver(x,y,Vx*ScaleV,Vy*ScaleV,0,'r','linewidth',2);
hold off
axis([-2*AtomSpacing,(LAtoms+2)*AtomSpacing,...
    -2*AtomSpacing,(WAtoms+2)*AtomSpacing]);

subplot(2,2,2),plot(time,KETot,'linewidth',2);
hold on
subplot(2,2,2),plot(time,PhiTot+LJEpsilon*nAtoms,'g','linewidth',2);
subplot(2,2,2),plot(time,PhiTot+KETot+LJEpsilon*nAtoms,'k','linewidth',2);
%     axis([0 TStop 0e-22 3e-23]);

AE = 1/2*Mass0*V2 + Phi-Phi0;
if AddParticle, AE(end) = 0;end
% AE = Phi;
subplot(2,2,4),scatter3(x,y,AE,ones(1,nAtoms)*300,AE,'fill');
axis([-2*AtomSpacing,(LAtoms+2)*AtomSpacing,...
    -2*AtomSpacing,(WAtoms+2)*AtomSpacing]);
view(2);

end

