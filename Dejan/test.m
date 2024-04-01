
F1 = @(x,y) -24*(x - y).^2;
Gc1 = @(x) x.^4;         %zgoraj
Gd1 = @(x) (x-2).^4;     %spodaj
Ga1 = @(y) y.^4;         %levo
Gb1 = @(y) (y-1).^4;     %desno


%1.a
[U1a,X1a,Y1a,A1a,B1a]= resiPoissonDif(0,1,0,2,F1,Ga1,Gb1,Gc1,Gd1,1,3);




%1.b
[U1b,X1b,Y1b,A1b,B1b]= resiPoissonDif(0,1,0,2,F1,Ga1,Gb1,Gc1,Gd1,3,3);
%--------------------------------------------------------------------


F2 = @(x,y) 4.*ones(size(x));
Gc2 = @(x) -x.^2;         %zgoraj
Gd2 = @(x) -x.^2 - 1;     %spodaj
Ga2 = @(y) -y.^2;         %levo
Gb2 = @(y) -y.^2 - 1;     %desno

%2.a
[U2a,X2a,Y2a,A2a,B2a]= resiPoissonDif(0,1,0,1,F2,Ga2,Gb2,Gc2,Gd2,3,3);

%2.b
[U2b,X2b,Y2b,A2b,B2b]= resiPoissonDif(0,1,0,1,F2,Ga2,Gb2,Gc2,Gd2,7,7);

%Primerjava vrednosti v točkah (1/4 2/4 3/4) x (1/4 2/4 3/4)
P2a = U2a(2:end-1,2:end-1);
P2b = U2b(3:2:end-2,3:2:end-2);

P = abs(P2a - P2b); % Vrednosti v istih točkah so pri obeh razmikih enake.

%--------------------------------------------------------------------------
%3.

% točna rešitev
tocnares = @(x,y) (x-y).^4;

% Razmik, ki je v obeh smereh enak postopoma zmanjšujemo -> število
% notranjih točk se povečuje
max_razlika = [];
for i = 1:10
    %figure
    [U,X,Y,A,B]= resiPoissonDif(0,1,0,2,F1,Ga1,Gb1,Gc1,Gd1,i,2*i+1);
    T = tocnares(X,Y);
    razlika = abs(U - T);
    max_razlika(i) = max(max(razlika)); %max_razlika konvergira proti 0, ko manjšamo razmika
    surf(X,Y,razlika);
end

