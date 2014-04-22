% clear all
clearvars
clearvars -GLOBAL
close all
format shorte
set(0,'DefaultFigureWindowStyle','docked')
global C V Mun Mup Gv Dn Dp Bv Em DnM MunM DpM MupM np pp x xm n p
global Rho divFp divFn niSi TwoCarriers t EpiSi n0 p0 tauSi
global Coupled NetDoping l PlotYAxis PlotSS im map
global PlotCount doPlotImage

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb*2*pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.Mun_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;% speed of light

Temp = 300;
C.Vt = C.kb*Temp/C.q_0;

EpiSi = C.eps_0*11.68;
MunSi = 1400*1e-4; % cm2 V-1s-1 * 1/(100 cm/m)^2
DnSi = MunSi*C.kb*Temp/C.q_0; % D = kt/q Mun
MupSi = 450*1e-4; % cm2 V-1s-1 * 1/(100 cm/m)^2
DpSi = MunSi*C.kb*Temp/C.q_0; % D = kt/q Mun
tauSi = 1e-8;
niSi = 1e10*1e6; % 1/cm^3 * (100 cm/m)^3 intrinsic concentration

JBC = 0; % No flow BC by default
RVbc = 0; % Ground Rightside 
SecondSim = 0;

PlotSS = 1;
PlotFile = 'image.gif'
PlotCount = 0;
doPlotImage = 0;

Simulation = 'PNJctEqBias'

if strcmp(Simulation,'GaussianTwoCar')
    eval('SetGaussian2CarParas');
elseif strcmp(Simulation,'GaussianTwoCarRC')
    eval('SetGaussian2CarParasRCOnly');
elseif strcmp(Simulation,'GaussianSingle0V')
    eval('SetGaussian1CarParas0V');
elseif strcmp(Simulation,'GaussianSingle1V')
    eval('SetGaussian1CarParas1V');
elseif strcmp(Simulation,'ExpDoping')
    eval('SetExpDopingParas');
elseif strcmp(Simulation,'PNJct')
    eval('SetPNJctParas');
elseif strcmp(Simulation,'PNJctEq')
    eval('SetPNJctParasEqBC');
elseif strcmp(Simulation,'PNJctEqBias')
    eval('SetPNJctParasEqBCBias');
end

FormGv(nx,LVbc,RVbc); % Poisson equation set Gv and Bv
[L,U] = lu(Gv);

Mun = ones(1,nx)*MunSi;
Dn = ones(1,nx)*DnSi;
MunM(1:nx-1) = (Mun(1:nx-1) + Mun(2:nx))/2;
DnM(1:nx-1) = (Dn(1:nx-1) + Dn(2:nx))/2;
n = zeros(1,nx);

Mup = ones(1,nx)*MupSi;
Dp = ones(1,nx)*DnSi;
MupM(1:nx-1) = (Mup(1:nx-1) + Mup(2:nx))/2;
DpM(1:nx-1) = (Dp(1:nx-1) + Dp(2:nx))/2;
p = zeros(1,nx);

if TwoCarriers == 1
    ni = NetDoping >= 0;
    n0(ni) = (NetDoping(ni) + sqrt(NetDoping(ni).^2 + 4* niSi*niSi))/2;
    p0(ni)  = niSi^2./n0(ni);
    
    pi = ~ni;
    p0(pi) = (-NetDoping(pi) + sqrt(NetDoping(pi).^2 + 4* niSi*niSi))/2;
    n0(pi)  = niSi^2./p0(pi);        
else
    n0  = NetDoping;
    p0 = zeros(1,nx);
end

n0 = n0 + npDisturbance;
if TwoCarriers == 1
    p0 = p0 + npDisturbance;
end

divFn = zeros(1,nx);
divFp = zeros(1,nx);

Rho = zeros(1,nx);
if (Coupled)
    Rho = C.q_0*(NetDoping - n0 + p0); % update Rho
    
    Rho(1) = 0; % 1 and nx are BC's
    Rho(nx) = 0;
end

V = U\(L\(-dx^2/EpiSi*Rho' + Bv'));
Em(1:nx-1) = -(V(2:nx) - V(1:nx-1))/dx;
MaxEm = max(abs(Em));
Maxn = max(n0);
Ld = sqrt(EpiSi/(C.q_0*Maxn));
dxMax = Ld/5;

dtMax = min(dx^2/2/max(Dn),dx^2/2/max(Dp));

if MaxEm > 0
    dt = min([2*dx/MaxEm dtMax])/4;
else
    dt = dtMax/4
end

t = 0;

n = n0;
np = n0;

p = p0;
pp = p0;

PlotVals(nx,dx,'on',l,TStop,PlotYAxis);

SimulateFlow(TStop,nx,dx,dtMax,JBC,RC,U,L,PlDelt)

if SecondSim == 1
    
   FormGv(nx,LVbc2,RVbc); % Poisson equation set Gv and Bv
   [L,U] = lu(Gv);

    SimulateFlow(TStop2,nx,dx,dtMax,JBC,RC,U,L,PlDelt)
end


if doPlotImage
    imwrite(im,map,PlotFile,'DelayTime',0.2,'LoopCount',inf);
end


fig2 = figure('Position', [100, 100, 1049, 895]);
PlotVals(nx,dx,'off',l,TStop,[]);

if doPlotImage
    f = getframe(fig2);
    [im,map] = rgb2ind(f.cdata,256,'nodither');

    filename = strcat('final-',PlotFile);
    imwrite(im,map,filename);
end
