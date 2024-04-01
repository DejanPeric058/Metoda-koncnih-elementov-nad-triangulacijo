p =@(x,y) 0*x + 1;
q =@(x,y) 0*x + 1;
r =@ (x,y) 0*x;
f =@(x,y) 0*x + 1;
g = @(x,y) fun_g(x,y);


[X, Y] = meshgrid(linspace(0, 1, 11));
X = X(:);
Y = Y(:);
TRI = delaunay(X, Y);
t = triangulation(TRI, X, Y);
tP = t.Points;
tC = t.ConnectivityList;

res = mke(p,q,r,f,t,g);
resP = res.Points;
resC = res.ConnectivityList;
trisurf(resC,resP(:,1), resP(:,2),resP(:,3))