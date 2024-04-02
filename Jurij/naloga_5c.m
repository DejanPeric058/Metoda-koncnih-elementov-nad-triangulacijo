%%
L = load('slo.mat');
t = triangulation(L.TRI, L.X, L.Y);
tP = t.Points;
tC = t.ConnectivityList;
triplot(t) 

p = @(x,y) -x;
q = @(x,y) -y;
r = @(x,y) 4*pi^2*(x + y);
f = @(x,y) 2*pi*cos(2*pi*(x + y));
g = @(x,y) cos(2*pi*x).*sin(2*pi*y);

res = mke(t,p,q,r,f,g);
resP1 = res.Points;
resC1 = res.ConnectivityList;

trisurf(resC1,resP1(:,1), resP1(:,2),resP1(:,3)); 
%%
res_actual = arrayfun(g, tP(:, 1), tP(:, 2));
trisurf(resC1, resP1(:,1), resP1(:,2), res_actual)
%%
trisurf(resC1, resP1(:,1), resP1(:,2), abs(resP1(:,3)-res_actual))
error = max(abs(resP1(:,3)-res_actual))
