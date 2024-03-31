function u = mke(p,q,r,f,t,g)
% Opis:
%  mke izracuna priblizek za resitev parcialne diferencialne enacbe
%   - d/dx (p(x,y) du/dx) - d/dy (q(x,y) du/dy) + r(x,y) u = f(x,y)
%  z robnim pogojem u = g po metodi koncnih elementov z zveznimi odsekoma
%  linearnimi funkcijami nad triangulacijo
%
% Definicija
%  u = mke(p,q,r,f,t,g)
%
% Vhodni podatki:
%  p,q,r,f  funkcije dveh spremenljivk, ki dolocajo parcialno diferencialno
%           enacbo,
%  t        triangulacija obmocja, predstavljena z razredom triangulation,
%  g        funkcija dveh spremenljivk, ki doloca vrednost resitve na robu
%           obmocja
%
% Izhodni podatek:
%  u        3D triangulacija, predstavljena z razredom triangulation, ki
%           doloca zvezno odsekoma linearno funkcijo, ki je priblizek za
%           resitev robnega problema po metodi najmanjsih kvadratov
P = t.Points;
T = t.ConnectivityList;
n = length(P);

F = freeBoundary(t);

neRobni = []; % točke, ki niso na robu
Robni = []; % točke, na robu 
for i = 1:n
    if ~ismember(i,F)
        neRobni = [neRobni ,i];
    else
        Robni = [Robni,i];
    end
end

ploscina =@(P1,P2,P3) 1/2*abs(P1(1)*(P2(2) - P3(2)) + P2(1)*(P3(2) - P1(2)) + P3(1)*(P1(2)-P2(2)));
% ta zanka določi b
m = length(neRobni);
b = zeros(m,1);
for i = 1:m
    bi = 0;
    for j = 1:length(T) % grem po vseh trikotnikih
        [bool,ID] =ismember(neRobni(i),T(j,:));
        if bool
            visine = zeros(3,1);
            visine(ID) = 1;
            Ptrikotnik = zeros(3,2);
            for tocka = 1:3
                Ptrikotnik(tocka,:) = P(T(j,tocka),:);
            end
            integral = triintegral(@(x,y)f(x,y).*trilin(Ptrikotnik,visine,x,y,'o'),Ptrikotnik);
            bi = bi + integral;
        end
    end
    % naslednja koda izracuna del integralov ki vsebuje odvode po x in y 
    sosedi = [];
    for trikotnik= T'
        if ismember(neRobni(i),trikotnik)
            sosedi = [sosedi, trikotnik];
        end
    end
    for trikotnik = sosedi
        for j = trikotnik'
            [bool,~] =ismember(j,Robni);
            if bool
                % naslednje vrstie so pomembne če je trikotnik na robu
                koordTrik = [P(trikotnik(1),:); P(trikotnik(2),:); P(trikotnik(3),:)];
                S = ploscina(koordTrik(1,:), koordTrik(2,:), koordTrik(3,:));
                
                % kot v notranjosti
                [~,ID] =ismember(neRobni(i),trikotnik);
                visina_i = zeros(3,1);
                visina_i(ID) = 1; 
        
                % kota na robu
                [~,JD] =ismember(j,trikotnik);
                visina_j = zeros(3,1);
                visina_j(JD) = g(P(j,1),P(j,2));
            
                po_1x = trilin(koordTrik,visina_i,1,1,'x');
                po_1y = trilin(koordTrik,visina_i,1,1,'y');
                po_2x = trilin(koordTrik,visina_j,1,1,'x');
                po_2y = trilin(koordTrik,visina_j,1,1,'y');
                
                int_p = triintegral(@(x,y)p(x,y),koordTrik);
                int_q = triintegral(@(x,y)q(x,y),koordTrik);
            
                bi = bi + S*(int_p*po_1x*po_2x + int_q*po_1y*po_2y);
            end 
        end
    end
    b(i) = bi;
end

% ta zanka je za določit A
A = zeros(m); % m je stevilo notranjih točk
for i = 1:m
    sosedi = [];
    for trikotnik= T'
        if ismember(neRobni(i),trikotnik)
            sosedi = [sosedi, trikotnik];
        end
    end
    for trikotnik = sosedi
        koordTrik = [P(trikotnik(1),:); P(trikotnik(2),:); P(trikotnik(3),:)];
        S = ploscina(koordTrik(1,:), koordTrik(2,:), koordTrik(3,:));
        
        [~,ID] =ismember(neRobni(i),trikotnik);
        visina_i = zeros(3,1);
        visina_i(ID) = 1; 
        for j = trikotnik'
            [bool,IR] =ismember(j,neRobni);
            if bool
                [~,JD] =ismember(j,trikotnik);
                visina_j = zeros(3,1);
                visina_j(JD) = 1;
            
                po_1x = trilin(koordTrik,visina_i,1,1,'x');
                po_1y = trilin(koordTrik,visina_i,1,1,'y');
                po_2x = trilin(koordTrik,visina_j,1,1,'x');
                po_2y = trilin(koordTrik,visina_j,1,1,'y');
                
                int_p = triintegral(@(x,y)p(x,y),koordTrik);
                int_q = triintegral(@(x,y)q(x,y),koordTrik);
            
                A(i,IR) = A(i,IR) + S*(int_p*po_1x*po_2x + int_q*po_1y*po_2y);
            end 
        end
    end
end
%min(min(A' == A))
alfe =  A\b;
P_res = [P,zeros(n,1)];
for i = 1:m
    P_res(neRobni(i),3) = alfe(i);
end
for i = 1:length(Robni)
    P_res(Robni(i),3) = g(P_res(Robni(i),1), P_res(Robni(i),2))
end

u = triangulation(T,P_res);
end




