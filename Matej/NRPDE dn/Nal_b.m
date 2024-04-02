p = @(x,y) 0*x + 1;
q = @(x,y) 0*x + 1;
r = @(x,y) 0*x + 0;
f = @(x,y) 0*x + 1;
%syms x y
%g = piecewise((0<=x)&(x<=1)&(y==0), x^3,(0<=x)&(x<=1)&(y==1),x^2,(x==0)&(0<=y)&(y<=1),sin(2*pi*y),(x==0)&(0<=y)&(y<=1),cos(2*pi*y));
g = @(x,y)funk_g(x,y);

[X,Y] = meshgrid(linspace(0,1,11));
X = X(:); Y = Y(:);
TRI = delaunay(X,Y);
T = triangulation(TRI,X,Y);

u = mke_vaje(p,q,r,f,g,T)
%class(g)
