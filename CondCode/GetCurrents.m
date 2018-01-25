function [ Curr ] = GetCurrents(ngeo, Max, nx, ny, Acond, Bcond, plot, SimType,FigNum)
global im fig fc map C

if SimType == 'c'
    pgeo = rand(ngeo, 2);
    pgeo(:, 1) = pgeo(:, 1) * nx;
    pgeo(:, 2) = pgeo(:, 2) * ny;
    rgeo = rand(ngeo,1)*Max;
    rcircs2 = rgeo.^2;

    cMap = zeros(ny,nx);

    for i = 1:nx
        for j = 1:ny
            cMap(j,i) = Acond;
            for p = 1:ngeo
                dx = (pgeo(p, 1) - i);
                dy = (pgeo(p, 2) - j);
                if (dx^2 + dy^2) < rcircs2(p)
                    cMap(j, i) = Bcond;
                end
            end
        end
    end
elseif SimType == 'r'
    pgeo = rand(ngeo,2);
    pgeo(:,1) = pgeo(:,1)*nx;
    pgeo(:,2) = pgeo(:,2)*ny;
    dgeo = ones(ngeo,1)*Max;
    
    cMap = zeros(ny,nx);
    
    for i = 1:nx
        for j = 1:ny
            cMap(j,i) = Acond;
            for p = 1:ngeo
                dx = abs(pgeo(p, 1) - i);
                dy = abs(pgeo(p, 2) - j);
                if (dx < dgeo(p)) && (dy < dgeo(p))
                    cMap(j, i) = Bcond;
                end
            end
        end
    end
else %if SimType == 'e'
    pgeo = rand(ngeo, 3);
    pgeo(:, 1) = pgeo(:, 1) * nx;
    pgeo(:, 2) = pgeo(:, 2) * ny;
    pgeo(:, 3) = pgeo(:, 2) * 180;
    rgeo = ones(ngeo, 1) * Max;

    cMap = zeros(ny, nx);

    for i = 1:nx
        for j = 1:ny
            cMap(j,i) = Acond;
            for p = 1:ngeo
                point = [j,i];
                el = [pgeo(p, 1), pgeo(p, 2),rgeo(p), rgeo(p) / 4, pgeo(p, 3)];
                if isPointInEllipse(point, el)
                    cMap(j,i)= Bcond;
                end
            end
        end
    end
end

G = sparse(nx*ny);
B = zeros(1,nx*ny);

for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        if i == 1
            G(n, :) = 0;
            G(n, n) = 1;
            B(n) = 1;
        elseif i == nx
            G(n, :) = 0;
            G(n, n) = 1;
        elseif j == 1
            nxm = j + (i - 2) * ny;
            nxp = j + (i) * ny;
            nyp = j + 1 + (i - 1) * ny;

            rxm = (cMap(j,i) + cMap(j,i-1)) / 2.0;
            rxp = (cMap(j,i) + cMap(j,i+1)) / 2.0;
            ryp = (cMap(j,i) + cMap(j+1,i)) / 2.0;

            G(n, n) = -(rxm+rxp+ryp);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nyp) = ryp;

        elseif j ==  ny
            nxm = j + (i - 2) * ny;
            nxp = j + (i) * ny;
            nym = j - 1 + (i - 1) * ny;

            rxm = (cMap(j,i) + cMap(j,i-1)) / 2.0;
            rxp = (cMap(j,i) + cMap(j,i+1)) / 2.0;
            rym = (cMap(j,i) + cMap(j-1,i)) / 2.0;

            G(n, n) = -(rxm + rxp + rym);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(j,i) + cMap(j,i-1))/2.0;
            rxp = (cMap(j,i) + cMap(j,i+1))/2.0;
            rym = (cMap(j,i) + cMap(j-1,i))/2.0;
            ryp = (cMap(j,i) + cMap(j+1,i))/2.0;

            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
        end

    end
end

V = G\B';

Vmap = reshape(V,[ny,nx]);

[Ex,Ey] = gradient(Vmap,1,1);

Ex = -Ex;
Ey = -Ey;

eFlowx = cMap .* Ex;
eFlowy = cMap .* Ey;

if plot
    if fc == 1
        fig = figure();
    end
    figure(FigNum)
    subplot(2, 2, 1), H = surf(cMap);
    set(H, 'linestyle', 'none');
    view(0, 90)
    subplot(2, 2, 2), H = surf(Vmap);
    set(H, 'linestyle', 'none');
    view(0, 90)
    subplot(2, 2, 3), quiver(Ex, Ey);
    axis([0 nx 0 ny]);
    subplot(2, 2, 4), quiver(eFlowx, eFlowy);
    axis([0 nx 0 ny]);
    
    if fc == 1
        f = getframe(fig);
        [im, map] = rgb2ind(f.cdata, 256, 'nodither');
        im(1, 1, 1, 2) = 0;
    else
        f = getframe(fig);
        im(:, :, 1, fc) = rgb2ind(f.cdata, map, 'nodither');
    end
end

C0 = sum(eFlowx(:,1));
Cnx = sum(eFlowx(:,ny));
Curr = (C0 + Cnx) * 0.5;

end