%%
p = @(x,y) 0*x + 1;
q = @(x,y) 0*x + 1;
r = @ (x,y) 0*x;
f = @(x,y) 0*x + 1;
g = @(x,y) fun_g(x,y);

[X, Y] = meshgrid(linspace(0, 1, 11));
X = X(:);
Y = Y(:);
TRI = delaunay(X, Y);
t = triangulation(TRI, X, Y);
tP = t.Points;
tC = t.ConnectivityList;

res = mke(t,p,q,r,f,g);
resP = res.Points;
resC = res.ConnectivityList;
trisurf(resC,resP(:,1), resP(:,2),resP(:,3))
colorbar
%%
fl = @(y) sin(2.*pi.*y);
fr = @(y) cos(2.*pi.*y);
fd = @(x) x.^3;
fu = @(x) x.^2;
f = @(x,y) x.*0 + 1;
Z2 = resiPoissonDif(0,1,0,1,f,fl, fr, fd, fu, 9, 9);
Z2 = reshape(Z2, 1, []);
error = max(abs(resP(:,3)'-Z2))