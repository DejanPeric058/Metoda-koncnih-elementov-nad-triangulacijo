
T = [1,2,9;1,9,8;8,9,7;9,6,7;9,10,6;2,9,10;2,3,10;3,4,10;10,4,5;10,5,6];    % Ogljisca trikotnikov
P = [0,0;1/2,0;1,0;1,1/2;1,1;1/2,1;0,1;0,1/2;1/4,1/2;3/4,1/2];              % Koordinate Tock

TR = triangulation(T,P);            % Definiramo triangulacijo
TR.ConnectivityList;                % Seznam povezav
TR.Points;                          % Seznam tock

triplot(T,P(:,1),P(:,2))            % Po potrebi narisemo triangulacijo

T = [0 0; 1 0; 0 1];                % Tocke
t = [0; 2; -1];                     % Vrednosti
[X,Y] = meshgrid(0:0.25:1,0:0.25:1);    

trilin(T,t,X,Y,'o');                % Poracunamo vrednosti (testiramo)

triintegral(@(x,y)trilin(T,t,x,y,'y'),T);       % Resimo integral po obmocju (testiramo)

f = @(x,y)cos(x.^2).*sin(3*y);      % Definiramo funkcijo za naprej
triintegral(f,T);                   % Poskusimo integrirati

% Zacetni in robni pogoji 

p =@(x,y) 0*x + 1;
q =@(x,y) 0*x + 1;
r =@ (x,y) 0*x;
f =@(x,y) 0*x + 1;
g = @(x,y) 0;

% res = mke_vaje(p,q,r,f,TR,g);
% resP = res.Points;
% resC = res.ConnectivityList;
% trisurf(resC,resP(:,1), resP(:,2),resP(:,3))

% Naredimo triangulacijo
[X, Y] = meshgrid(linspace(0, 1, 11));
X = X(:);
Y = Y(:);
TRI = delaunay(X, Y);
t = triangulation(TRI, X, Y);
tP = t.Points;
tC = t.ConnectivityList;
% triplot(tC,tP(:,1), tP(:,2))


% res = mke_vaje(p,q,r,f,t,g);
% resP = res.Points
% resC = res.ConnectivityList;
%trisurf(resC,resP(:,1), resP(:,2),resP(:,3))

g = @(x,y)fun_g(x,y);                % Robni pogoji
% KLICEMO FUNKCIJO
res = mke_vaje(p,q,r,f,t,g);
resP1 = res.Points;
resC1 = res.ConnectivityList;
trisurf(resC1,resP1(:,1), resP1(:,2),resP1(:,3));       % Samo narisemo
%%%% REZULTATE BEREMO KOT ZADNJI STOLPEC V RES.POINTS 

% Naredimo ven matriko
K = 11;
J = 11;
matriziraj_vec =@(vec) rot90(flip(reshape(flip(vec),J,K)',2),3);
mat_KoncnihElementov = matriziraj_vec(resP1(:,3));

% Primerjamo z diferencami
f =@(x,y) 0*x + 1;

J = 19;
K = 19;
a = 0;
b = 1;
c = 0;
d = 1;

Gc = @(x) x.^3;
Gd = @(x) x.^2;
Ga = @(y) sin(2*pi*y);
Gb = @(y) cos(2*pi*y);

[U,X,Y,A,B] = resiPoissonDif(a,b,c,d,f,Gc,Gd,Ga,Gb,J,K);
diference = flip(U);

max(max(abs(diference - mat_KoncnihElementov)));


% Nalozimo triangulacijo iz predloge
L = load('slo.mat');
t = triangulation(L.TRI, L.X, L.Y);

% pazi na minuse, v originalni so minusi
p =@(x,y) -x;
q =@(x,y) -y;
r =@ (x,y) 4*pi^2*(x + y);
f =@(x,y) 2*pi*cos(2*pi*(x + y));
g = @(x,y) cos(2*pi*x).*sin(2*pi*y);

res = mke_vaje(p,q,r,f,t,g);
resP1 = res.Points;
resC1 = res.ConnectivityList;
% Izpustimo, lahko preberemo vse iz resP1
% mat_KoncnihElementov = matriziraj_vec(resP1(:,3));



u = @(x,y) cos(2*pi*x).*sin(2*pi*y);
max_napaka = 0;
for P = resP1'
    napaka = norm(P(3) - u(P(1),P(2)));
    if max_napaka < napaka
        max_napaka = napaka;
    end
end
max_napaka



