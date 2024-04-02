function [f, f_dx, f_dy] = trilin(T,t)
% Opis:
%  trilin izracuna vrednosti (odvodov) linearne funkcije, ki je dolocena z
%  vrednostmi v ogliscih trikotnika
%
% Definicija:
%  [f, f_dx, f_dy] = trilin(T,t)
%
% Vhodni podatki:
%  T    tabela velikosti 3 x 2, v kateri vsaka vrstica predstavlja
%       kartezicni koordinati oglisca trikotnika
%  t    stolpec dolzine 3, v katerem vsak element predstavlja vrednost
%       linearne funkcije v ogliscu trikotnika,
% Izhodni podatki:
%  f_dx    vrednosti odvoda linearne funkcije po x
%  f_dy    vrednosti odvoda linearne funkcije po y
%  f       predpis za linearno funckijo

sistem = [ones(3,1) T];
koef = sistem \ t;

f_dx = koef(2);
f_dy = koef(3);
f = @(x,y) koef(1) + koef(2)*x + koef(3)*y;
end


