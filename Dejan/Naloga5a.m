
p =@(x,y) 0*x + 1;
q =@(x,y) 0*x + 1;
r =@ (x,y) 0*x;
f =@(x,y) 0*x + 1;
g = @(x,y) 0;

T = [1,4,5;1,2,5;2,5,6;2,3,6;3,6,7;4,5,8;5,8,9;5,6,9;6,9,10;6,7,10];    % Ogljisca trikotnikov
X = [0; 0.5; 1; 0; 0.25; 0.75; 1; 0; 0.5; 1];
Y = [0; 0; 0; 0.5; 0.5; 0.5; 0.5; 1; 1; 1];


TR = triangulation(T,X,Y); 

res = mke(p,q,r,f,TR,g);
resP = res.Points;
resC = res.ConnectivityList;
trisurf(resC,resP(:,1), resP(:,2),resP(:,3))
colormap([0 0 1])