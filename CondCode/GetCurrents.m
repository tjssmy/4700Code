function [ Curr ] = GetCurrents(ncircs,MaxRad,nx,ny,Acond,Bcond,plot)

pcircs = rand(ncircs,2);
pcircs(:,1) = pcircs(:,1)*nx;
pcircs(:,2) = pcircs(:,2)*ny;
rcircs = rand(ncircs,1)*MaxRad;
rcircs2 = rcircs.^2;

cMap = zeros(nx,ny);

for i = 1:nx
    for j = 1:ny
        cMap(i,j) = Acond;
        for p = 1:ncircs
            dx = (pcircs(p,1)-i);
            dy = (pcircs(p,2)-j);
            if dx*dx + dy*dy < rcircs2(p)
                cMap(i,j)= Bcond;
            end
        end
    end
end


G = sparse(nx*ny);
B = zeros(1,nx*ny);

for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        
        if i == 1 
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = 1;
        elseif i == nx
            G(n,:) = 0;
            G(n,n) = 1;
        elseif j == 1
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nyp = j+1 + (i-1)*ny;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;
            
            G(n,n) = -(rxm+rxp+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;
            
        elseif j ==  ny
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
                        
            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
                        
            G(n,n) = -(rxm+rxp+rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
                        
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;
            
            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;
            
            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
        end
        
    end
end



V = G\B';

Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        
        Vmap(i,j) = V(n);
    end
end

for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i,j) = (Vmap(i+1,j) - Vmap(i,j));
        elseif i == nx
            Ex(i,j) = (Vmap(i,j) - Vmap(i-1,j));
        else 
            Ex(i,j) = (Vmap(i+1,j) - Vmap(i-1,j))*0.5;
        end
        if j == 1
            Ey(i,j) = (Vmap(i,j+1) - Vmap(i,j));
        elseif j == ny
            Ey(i,j) = (Vmap(i,j) - Vmap(i,j-1));
        else 
            Ey(i,j) = (Vmap(i,j+1) - Vmap(i,j-1))*0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;

eFlowx = cMap.*Ex;
eFlowy = cMap.*Ey;

if plot
    subplot(2,2,1),surface(cMap');
    subplot(2,2,2),surface(Vmap');
    subplot(2,2,3),quiver(Ex',Ey');
    axis([0 nx 0 ny]);
    subplot(2,2,4),quiver(eFlowx',eFlowy');
    axis([0 nx 0 ny]);
end

C0 = sum(eFlowx(1,:));
Cnx = sum(eFlowx(nx,:));

Curr = (C0 + Cnx)*0.5;


end

