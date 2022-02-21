function SimulateFlow(TStop,nx,dx,dtMax,JBC,RC,U,L,PlDelt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global C V  Bv Em np pp  n p n0 p0 Coupled
global Rho divFp divFn niSi TwoCarriers t EpiSi tauSi
global NetDoping l  PlotYAxis Urc

Plt0 = PlDelt;

t0 = t;
while t < TStop

    MaxEm = max(abs(Em));

    if MaxEm > 0
        dt = min([2*dx/MaxEm dtMax])/4;
    else
        dt = dtMax/4;
    end

    V = U\(L\(-dx^2/EpiSi*Rho' + Bv'));
    Em(1:nx-1) = -(V(2:nx) - V(1:nx-1))/dx;

    DivFn(nx,dx,JBC);
    n = np + dt*divFn;

    Maxn = max(n);

    if TwoCarriers == 1
        DivFp(nx,dx,JBC);
        p = pp + dt*divFp;
        Maxp = max(p);
        if RC
            Ecrit = 5e6;
            iE = abs(Em) > Ecrit;
           
            tauSiv = tauSi*ones(size(n));
%             tauSiv(iE) = 1e-6;
            Del = zeros(size(n));
            Del(iE) = (n(iE)+p(iE)).*exp(Em(iE)/Ecrit)*1e-5;
            Urc = (n.*p - n0.*p0)./(n + p + 2*niSi)./tauSiv;
            dn = zeros(size(n));
            if t > 5.28001e-9 && t  < 5.28002e-9                
                dn(iE) = 1e20;
            end
                
            n = n - dt*Urc+Del+dn;
            p = p - dt*Urc+Del+dn;
        end
    else
        Maxp = 0;
    end

    Ldn = sqrt(EpiSi/(C.q_0*Maxn));

    if TwoCarriers == 1
        Ldp = sqrt(EpiSi/(C.q_0*Maxp));
    else
        Ldp = 1e6;
    end

    dxMax = min(Ldn,Ldp)/5;

    if (Coupled)
        Rho = C.q_0*(NetDoping - n + p); % update Rho
        Rho(1) = 0; % 1 and nx are BC's
        Rho(nx) = 0;
    end


    np = n;
    pp = p;

    t = t + dt;
    if t > Plt0
        fprintf('time: %g (%5.2g %%)\n',t,t/TStop*100);
        PlotValsSimple(nx,dx,'off',l,TStop,PlotYAxis,1);
        Plt0 = Plt0 + PlDelt;
        pause(0.0001)
    end
end

end
