
load('Rev4p284.mat')


% PlotValsSimple(nx,dx,'off',l,TStop,[],1);
% % plot(Phi-V(1),log(abs(mean(Itot(:,end)))),'*');hold on
% plot(Phi-V(1),mean(Itot(:,end)),'*');hold on
% 
FormGv(nx,LVbc2,RVbc); % Poisson equation set Gv and Bv
[L,U] = lu(Gv);
TStop2 = TStop2 +  8000000*1e-18;
Itot = [];
tv = [];
V1 = [];


SimulateFlow(TStop2,nx,dx,dtMax,JBC,RC,U,L,PlDelt)
PlotValsSimple(nx,dx,'off',l,TStop,[],1);
subplot(3,3,9);
%     plot(Phi-V(1),log(abs(mean(Itot(:,end)))),'*');hold on
plot(Phi-V(1),-mean(Itot(:,end)),'*');hold on

% fig2 = figure('Position', [100, 100, 1049, 895]);
PlotValsSimple(nx,dx,'off',l,TStop,[],1);
