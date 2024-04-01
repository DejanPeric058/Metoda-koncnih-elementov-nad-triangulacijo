
p =@(x,y) -x;
q =@(x,y) -y;
r =@ (x,y) 4*pi^2*(x+y);
f =@(x,y) 2*pi*cos(2*pi*(x+y));
g = @(x,y) cos(2*pi*x)*sin(2*pi*y);


L = load('slo.mat');
t = triangulation(L.TRI, L.X, L.Y);
tP = t.Points;
tC = t.ConnectivityList;

res = mke(p,q,r,f,t,g);
resP = res.Points;
resC = res.ConnectivityList;
trisurf(resC,resP(:,1), resP(:,2),resP(:,3))

