function dist = eu_dist(x, y)

xshape = size(x);
Nx = xshape(1);
Ndim = xshape(2);

yshape = size(y);
Ny = yshape(1);

dist  = zeros(Nx,Ny,Ndim);

for i = 1:Nx
    for j = 1:Ny
        for k = 1:Ndim
            dist(i,j,k) = x(i,k) - y(j,k);
        end
    end
end