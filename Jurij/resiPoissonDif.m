function [U,X,Y,A,B] = resiPoissonDif(a,b,c,d,F,Ga,Gb,Gc,Gd,J,K)
% Metoda za numerično reševanje Poissonove enačbe -laplace U = F na
% pravokotnem območju (a,b) x (c,d). 
%
% Vhodni podatki:
% a,b,c,d           parametri, ki določajo pravokotnik (a,b) x (c,d)
% F                 funkcija v enačbi -laplace U = F
% Ga, Gb, Gc, Gd    funkcije, ki določajo Dirichletove robne pogoje, kjer 
% Gc, Gd določata robne pogoje v x smeri:
%                   U(x,c) = Gc(x,c) 
%                   U(x,d) = Gd(x,d); za vsak x iz [a,b]
% Ga, Gb določata robne pogoje v y smeri:
%                   U(a,y) = Ga(a,y) 
%                   U(b,y) = Gb(b,y); za vsak y iz [c,d]
% J, K              parametra diskretizacije, ki določata število notranjih
%                   točk v x oziroma y smeri
%
% Izhodni podatki:
%  U        tabela velikosti (K+2) x (J+2), ki predstavlja numerično
%           rešitev Poissonove enacbe - laplace U = F pri danih robnih
%           pogojih,
%  X,Y      tabeli velikosti (K+2) x (J+2), ki vsebujeta x in y koordinate
%           točk mreže, U(j,k) predstavlja numerični približek za rešitev
%           Poissonove enačbe v točki (X(j,k), Y(j,k)),
%  A,B      matrika in vektor sistema A*u = B, katerega rešitev določa
%           numerične približke v notranjih točkah 

U = zeros(J+2,K+2); % Bomo kasneje transponirali, da bo (K+2) x (J+2)
% Diskretizacija pravokotnika:
deltax = (b-a)/(J+1);
deltay = (d-c)/(K+1);

x = [a:deltax:b];
y = [c:deltay:d];
[X,Y] = meshgrid(x,y);

delta_kvadrat = (deltax^2*deltay^2)/(2*(deltax^2 + deltay^2));
thetax = delta_kvadrat/deltax^2;
thetay = delta_kvadrat/deltay^2; % Velja thetax + thetay = 1/2

% Matrika A:
diag = ones(J*K,1);
poddiag1 = [-thetax*ones(J*K-1,1);0];
naddiag1 = [0;-thetax*ones(J*K-1,1)];
if J > 1
poddiag1(J:J:end) = 0;
naddiag1(J+1:J:end) = 0;
end
poddiag2 = [-thetay*ones(J*K-J,1);zeros(J,1)];
naddiag2 = [zeros(J,1);-thetay*ones(J*K-J,1)];

% spidiags ustvari razpršeno matriko velikost JK x JK iz podanih vektorjev,
% ki jih vstavi na diagonale podane v 2. argumentu funkcije. 
if J > 1
A = spdiags([poddiag2,poddiag1,diag,naddiag1,naddiag2],[-J,-1,0,1,J],J*K,J*K);
else
    A = spdiags([poddiag2,diag,naddiag2],[-1,0,1],J*K,J*K);
end

% Določimo še desno stran sistema. B = delta_kvadrat*f + g
f = F(X,Y);
f = f(2:end-1,2:end-1)';
B = delta_kvadrat*f(:);

B(1:J) = B(1:J) + thetay*Gc(X(1,2:end-1))';
B(end-J+1:end) =  B(end-J+1:end) + thetay*Gd(X(1,2:end-1))';
B(1:J:end) = B(1:J:end) + thetax*Ga(Y(2:end-1,1));
B(J:J:end) = B(J:J:end) + thetax*Gb(Y(2:end-1,1));

% Rešimo sistem enačb:
resitev = A\B;


U1 = reshape(resitev,[J,K]); % Vrne JK matriko z elementi rešitve
% Rešitve v notranjih točkah:
U(2:end-1,2:end-1) = U1;  
U = U'; % (K+2) x (J+2) matrika

% Dirichletovi robni pogoji:
U(1,:) = Gc(X(1,:));    %zgoraj
U(end,:) = Gd(X(1,:));  %spodaj
U(:,1) = Ga(Y(:,1));    %levo
U(:,end) = Gb(Y(:,1));  %desno

% Tridimenzionalen prikaz rešitve:
surf(X,Y,U);
end