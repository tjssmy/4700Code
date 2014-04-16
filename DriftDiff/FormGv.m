function  FormGv(nx,lBC)
global Gv Bv

Gv = sparse(nx,nx);
Bv = zeros(1,nx);

for i = 1:nx
    
    if i == 1
        if lBC == 'f'
            Gv(i,i) = 1;
            Gv(i,i-1) = -1; % floating contact
        else
            Bv(i) = lBC;
            Gv(i,i) = 1;
        end            
    elseif i == nx
        Gv(i,i) = 1;
    else
        Gv(i,i-1) = 1;
        Gv(i,i) = -2;
        Gv(i,i+1) = 1;
    end
end

end

