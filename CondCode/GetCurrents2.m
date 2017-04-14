function [ Curr ] = GetCurrents2(ngeo, Max, nx, ny, Acond, Bcond, plot, SimType)
global im fig fc map

ix = 1:nx;
iy = 1:ny;
[xx,yy] = meshgrid(ix,iy);
ft      = 0*xx;
cMap = zeros(nx,ny) + Acond;

if SimType == 'c'
    pgeo = rand(ngeo, 2);
    pgeo(:, 1) = pgeo(:, 1) * nx;
    pgeo(:, 2) = pgeo(:, 2) * ny;
    rgeo = rand(ngeo,1)*Max;
    rcircs2 = rgeo.^2;


    for p = 1:ngeo
        ft = ((pgeo(p,1) - xx).^2 + (pgeo(p,2) - yy).^2) <= rcircs2(p);
        cMap(ft) = Bcond;
    end

else %if SimType == 'e'
    pgeo = rand(ngeo, 3);
    pgeo(:, 1) = pgeo(:, 1) * nx;
    pgeo(:, 2) = pgeo(:, 2) * ny;
    pgeo(:, 3) = pgeo(:, 2) * 180;
    rgeo = ones(ngeo, 1) * Max;

    %pgeo(p,1) = xloc
    %pgeo(p,2) = yloc
    %rgeo(p,2) = xrad
    %rgeo(p,2) = yrad
    %pgeo(p,3) = rotation (degrees)
    % cMap = zeros(nx, ny);

    for p = 1:ngeo
        ft = (...
             (((xx-pgeo(p,1))*cosd(pgeo(p,3)) + (yy-pgeo(p,2))*sind(pgeo(p,3)))/(rgeo(p,1)*1)).^2 ...
           + (((xx-pgeo(p,1))*sind(pgeo(p,3)) - (yy-pgeo(p,2))*cosd(pgeo(p,3)))/(rgeo(p,1)/4)).^2 ...
             ) <= 1;
        cMap(ft) = Bcond;
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

            rxm = (cMap(i, j) + cMap(i - 1, j)) / 2.0;
            rxp = (cMap(i, j) + cMap(i + 1, j)) / 2.0;
            ryp = (cMap(i, j) + cMap(i, j + 1)) / 2.0;

            G(n, n) = -(rxm+rxp+ryp);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nyp) = ryp;

        elseif j ==  ny
            nxm = j + (i - 2) * ny;
            nxp = j + (i) * ny;
            nym = j - 1 + (i - 1) * ny;

            rxm = (cMap(i, j) + cMap(i - 1, j)) / 2.0;
            rxp = (cMap(i, j) + cMap(i + 1, j)) / 2.0;
            rym = (cMap(i, j) + cMap(i, j - 1)) / 2.0;

            G(n, n) = -(rxm + rxp + rym);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nym) = rym;
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
        n = j + (i - 1) * ny;

        Vmap(i, j) = V(n);
    end
end

for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i, j));
        elseif i == nx
            Ex(i, j) = (Vmap(i, j) - Vmap(i - 1, j));
        else
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i - 1, j)) * 0.5;
        end
        if j == 1
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j));
        elseif j == ny
            Ey(i, j) = (Vmap(i, j) - Vmap(i, j - 1));
        else
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j - 1)) * 0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;

eFlowx = cMap .* Ex;
eFlowy = cMap .* Ey;

if plot
    if fc == 1
        fig = figure();
    end

    subplot(2, 2, 1), H = surf(cMap');
    set(H, 'linestyle', 'none');
    view(0, 90); title('Current')
    subplot(2, 2, 2), H = surf(Vmap');
    set(H, 'linestyle', 'none');
    view(0, 90); title('Voltage')
    subplot(2, 2, 3), quiver(Ex', Ey');
    axis([0 nx 0 ny]); title('E Field')
    subplot(2, 2, 4), quiver(eFlowx', eFlowy');
    axis([0 nx 0 ny]); title('E field * Current')

    if fc == 1
        f = getframe(fig);
        [im, map] = rgb2ind(f.cdata, 256, 'nodither');
        im(1, 1, 1, 2) = 0;
    else
        f = getframe(fig);
        im(:, :, 1, fc) = rgb2ind(f.cdata, map, 'nodither');
    end
end

C0 = sum(eFlowx(1, :));
Cnx = sum(eFlowx(nx, :));
Curr = (C0 + Cnx) * 0.5;

end
