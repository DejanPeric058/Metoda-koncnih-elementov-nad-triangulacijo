function u = mke_vaje(p,q,r,f,g,t)
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
n = length(P);
C = t.ConnectivityList;

F = freeBoundary(t);

notranje = []; % točke, ki niso na robu
robne = []; % točke, na robu 

for i = 1:n
    if ~ismember(i,F)
        notranje = [notranje ,i];
    else
        robne = [robne,i];
    end
end

m = length(notranje);
A = zeros(m);
b = zeros(m,1);

for i = 1:size(C,1)
    
    T = C(i,:);
    trikotnik = [P(T(1),:); P(T(2),:); P(T(3),:)];
    
    for j = 1:3
        [je_notranja1,id1] = ismember(T(j),notranje);
        G1 = zeros(3,1);
        if je_notranja1
            G1(j) = 1;
        else
            G1(j) = g(P(T(j),1),P(T(j),2));
        end
            for k = j:3
                [je_notranja2,id2] = ismember(T(k),notranje);
                
                G2 = zeros(3,1);
                if je_notranja2
                    G2(k) = 1;
                else
                    G2(k) = g(P(T(k),1),P(T(k),2));
                end
           
                [p1,p1_x,p1_y] = trilin(trikotnik,G1);
                [p2,p2_x,p2_y] = trilin(trikotnik,G2);
           
                int_p = triintegral(@(x,y)p(x,y),trikotnik);
                int_q = triintegral(@(x,y)q(x,y),trikotnik);
                
                clen1 = int_p*p1_x*p2_x;
                clen2 = int_q*p1_y*p2_y;
                blabla = @(p) p*p1*p2;
                clen3 = triintegral(@(x,y) blabla(r(x,y)),trikotnik);
                
                if je_notranja2 + je_notranja1 == 2
                    A(id1,id2) = A(id1,id2) + clen1 + clen2 + clen3;
                    A(id2,id1) = A(id1,id2);
                    
                    if j == k
                        clen4 = triintegral(@(x,y) f(x,y).*p1,trikotnik);
                        b(id1) = b(id1) + clen4;
                    end
                else
                    if je_notranja1
                        b(id1) = b(id1) - clen1 - clen2 - clen3;
                    elseif je_notranja2
                        b(id2) = b(id2) - clen1 - clen2 - clen3; 
                    end
                end
            end
    end
end

alfe =  A\b;
P_u = [P,zeros(n,1)];

for i = 1:m
    P_u(notranje(i),3) = alfe(i);
end

for i = 1:length(robne)
    P_u(robne(i),3) = g(P_u(robne(i),1), P_u(robne(i),2));
end

u = triangulation(C,P_u);

end