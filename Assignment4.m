
%Part_D

for c=1:4     %alalytical approach
    dt=0.0002*10^(-c);
    sz=round(0.002/dt)+1;
    V=zeros(1,sz);
    T=zeros(1,sz);
    R=20;
    C=10e-6;
    V(1)=0;
    T(1)=0;
    for i=1:0.002/dt     %few data points (?)
        V(i+1)=(1+R*C*V(i)/dt)/(1+R*C/dt);
        T(i+1)=i*dt;
    end
    figure(c)
    plot (T,V)
end
