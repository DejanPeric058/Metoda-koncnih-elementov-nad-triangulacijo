function z = trilin(T,t,x,y,o)
% Opis:
%  trilin izracuna vrednosti (odvodov) linearne funkcije, ki je dolocena z
%  vrednostmi v ogliscih trikotnika
%
% Definicija:
%  z = trilin(T,t,x,y,o)
%
% Vhodni podatki:
%  T    tabela velikosti 3 x 2, v kateri vsaka vrstica predstavlja
%       kartezicni koordinati oglisca trikotnika
%  t    stolpec dolzine 3, v katerem vsak element predstavlja vrednost
%       linearne funkcije v ogliscu trikotnika,
%  x,y  seznama kartezicnih koordinat tock, v katerih racunamo vrednosti
%       funkcije,
%  o    parameter, ki doloca odvod: ce ni podan, metoda vraca vrednosti
%       linearne funkcije, ce je o = 'x' ali o = 'y' pa vrednosti odvodov
%       funkcije po x oziroma po y
%
% Izhodni podatek:
%  z    vrednosti (odvodov) linearne funkcije v tockah, dolocenih s
%       seznamoma x in y

sistem = [ones(3,1) T];

koef = sistem\t;

lin =@(x,y) koef(1) + koef(2)*x + koef(3)*y;

if o == 'x'
    z = koef(2)*(x*0 + 1);
elseif o == 'y'
    z = koef(3)*(y*0 + 1);
else
    z = lin(x,y);
end
end


