
load('Rev0p6.mat')


PlotValsSimple(nx,dx,'off',l,TStop,[],1);
% plot(Phi-V(1),log(abs(mean(Itot(:,end)))),'*');hold on
plot(Phi-V(1),mean(Itot(:,end)),'*');hold on

LVbc20 = LVbc2;

for k = 1:40
    LVbc2 = LVbc20+0.1*k;
%     tauSi = tauSi*10;
    FormGv(nx,LVbc2,RVbc); % Poisson equation set Gv and Bv
    [L,U] = lu(Gv);
    TStop2 = TStop2 +  80000000*1e-18;
    SimulateFlow(TStop2,nx,dx,dtMax,JBC,RC,U,L,PlDelt)
    PlotValsSimple(nx,dx,'off',l,TStop,[],1);
    subplot(3,3,9);
%     plot(Phi-V(1),log(abs(mean(Itot(:,end)))),'*');hold on
    plot(Phi-V(1),-mean(Itot(:,end)),'*');hold on
end

% fig2 = figure('Position', [100, 100, 1049, 895]);
PlotValsSimple(nx,dx,'off',l,TStop,[],1);
