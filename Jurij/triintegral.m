function integral = triintegral(f,T)
% Opis:
%  triintegral vrne priblizek za vrednost integrala funkcije na danem
%  trikotniku
%
% Definicija:
%  integral = triintegral(f,T)
%
% Vhodna podatka:
%  f    funkcija dveh spremenljivk,
%  T    tabela velikosti 3 x 2, v kateri vsaka vrstica predstavlja
%       kartezicni koordinati oglisca trikotnika, na katerem integriramo
%       funkcijo f
%
% Izhodni podatek:
%  integral     priblizek za integral funkcije f na trikotniku, podanem s tabelo T

P1 = T(1,:);
P2 = T(2,:);
P3 = T(3,:);

phi_x = @(u,v) u*P1(1) + v*P2(1) + (1-u-v)*P3(1);
phi_y = @(u,v) u*P1(2) + v*P2(2) + (1-u-v)*P3(2);
det_J_phi = (P1(1) - P3(1))*(P2(2) - P3(2)) - (P2(1) - P3(1))*(P1(2) - P3(2));

integral = integral2(@(u,v)f(phi_x(u, v),phi_y(u,v))*abs(det_J_phi),0,1,0,@(u)1-u);


