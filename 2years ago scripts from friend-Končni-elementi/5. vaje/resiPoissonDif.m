function [U,X,Y,A,B] = resiPoissonDif(a,b,c,d,F,Gc,Gd,Ga,Gb,J,K)
% Opis:
%  resiPoissonDif z uporabo diferencne metode resi Poissonovo enacbo na
%  pravokotniku pri Dirichletovih robnih pogojih
%
% Definicija:
%  [U,X,Y,A,B] = resiPoissonDif(a,b,c,d,F,Gc,Gd,Ga,Gb,J,K)
%
% Vhodni podatki:
%  a,b,c,d  parametri, ki dolocajo pravokotnik (a,b) x (c,d),
%  F        funkcija, ki doloca Poissonovo enacbo - laplace U = F,
%  Gc,Gd    funkciji, ki dolocata robne pogoje v x smeri:
%           U(x,c) = Gc(x), U(x,d) = Gd(x)
%  Ga,Gb    funkciji, ki dolocata robne pogoje v y smeri:
%           U(a,y) = Ga(y), U(b,y) = Gb(y)
%  J,K      parametra diskretizacije pri diferencni metodi, ki dolocata
%           stevilo notranjih tock mreze v x oziroma y smeri
%
% Izhodni podatki:
%  U        tabela velikosti (K+2) x (J+2), ki predstavlja numericno
%           resitev Poissonove enacbe - laplace U = F pri danih robnih
%           pogojih,
%  X,Y      tabeli velikosti (K+2) x (J+2), ki vsebujeta x in y koordinate
%           tock mreze, na kateri se izvede diferencna metoda (vrednost
%           U(k,j) torej predstavlja numericni priblizek za resitev
%           Poissonove enacbe v tocki (X(k,j), Y(k,j)),
%  A,B      matrika in vektor sistema A*x = B, katerega resitev doloca
%           numericne priblizke v notranjih tockah mreze

delta_x = (b - a)/(J+1);
delta_y = (d - c)/(K + 1);
delta = (delta_x^2 * delta_y^2)/(2*(delta_x^2 + delta_y^2));
theta_x = delta/delta_x^2;
theta_y = delta/delta_y^2;

C_diag = (1)*ones(1,J);
C_nadPodDiag = -theta_x*ones(1,J-1);
C = diag(C_diag) + diag(C_nadPodDiag,1) + diag(C_nadPodDiag,-1);

D_diag = -theta_y*ones(1,J);
D = diag(D_diag);

vzorcni = ones(1,K-1);
A = kron(eye(K),C) + kron(diag(vzorcni,1)+ diag(vzorcni,-1),D);

x = a:delta_x:b;
y = c:delta_y:d;

[X,Y] = meshgrid(x,y(end:-1:1));
[X_obrezan,Y_obrezan] = meshgrid(x(2:end-1),y(end-1:-1:2));
B_prviDel = F(X_obrezan,Y_obrezan);

mat_robnihY = [Ga(y(end-1:-1:2)'),zeros(K,J-2),Gb(y(end-1:-1:2)')];
mat_robnihX = [Gd(x(2:end-1));zeros(K-2,J);Gc(x(2:end-1))];

vektoriziraj_mat =@(A)reshape(flip(A)',K*J,1);

B_prviDel_vec = vektoriziraj_mat(B_prviDel);
B_drugiDel = vektoriziraj_mat(mat_robnihY);
B_tretjiDel = vektoriziraj_mat(mat_robnihX);

B = delta*B_prviDel_vec + theta_x*B_drugiDel + theta_y*B_tretjiDel;

U_notranji =A\B;

matriziraj_vec =@(vec) flip(reshape(flip(vec),J,K)',2);

U_notranjiMat = matriziraj_vec(U_notranji);

U_skoraj = [Ga(y(end-1:-1:2)'),U_notranjiMat,Gb(y(end-1:-1:2)')];
U = [Gd(x); U_skoraj; Gc(x)];
end















