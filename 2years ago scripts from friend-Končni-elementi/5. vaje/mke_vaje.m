function u = mke_vaje(p,q,r,f,t,g)
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

m = length(neRobni);
A = zeros(m);
b = zeros(m,1);

for k = 1:size(t.ConnectivityList,1)
    
    T = t.ConnectivityList(k,:);
    koordTrik = [P(T(1),:); P(T(2),:); P(T(3),:)];
    
    for I = 1:3
        [robI,IdI] = ismember(T(I),neRobni);
        visinaI = zeros(3,1);
        if robI
            visinaI(I) = 1;
        else
            visinaI(I) = g(P(T(I),1),P(T(I),2));
        end
            for J = I:3
                [robJ,IdJ] = ismember(T(J),neRobni);
                
                visinaJ = zeros(3,1);
                if robJ
                    visinaJ(J) = 1;
                else
                    visinaJ(J) = g(P(T(J),1),P(T(J),2));
                end
           
                po_1x = trilin(koordTrik,visinaI,1,1,'x');
                po_1y = trilin(koordTrik,visinaI,1,1,'y');
                po_2x = trilin(koordTrik,visinaJ,1,1,'x');
                po_2y = trilin(koordTrik,visinaJ,1,1,'y');
           
                int_p = triintegral(@(x,y)p(x,y),koordTrik);
                int_q = triintegral(@(x,y)q(x,y),koordTrik);
                
                prvi_del = int_p*po_1x*po_2x;
                drugi_del = int_q*po_1y*po_2y;
           
                rfiifij = @(x,y) r(x,y).*trilin(koordTrik,visinaJ,x,y,'o').*trilin(koordTrik,visinaI,x,y,'o');
                tretji_del = triintegral(rfiifij,koordTrik);
                if robJ + robI == 2
                    A(IdI,IdJ) = A(IdI,IdJ) + prvi_del + drugi_del + tretji_del;
                    A(IdJ,IdI) = A(IdI,IdJ);
                    
                    fi = @(x,y) trilin(koordTrik,visinaI,x,y,'o');
                    int_ffi = triintegral(@(x,y) f(x,y).*fi(x,y),koordTrik);
                    if I == J
                        b(IdI) = b(IdI) + int_ffi;
                    end
                else
                    if robI
                        b(IdI) = b(IdI) - prvi_del - drugi_del - tretji_del;
                    elseif robJ
                        b(IdJ) = b(IdJ) - prvi_del - drugi_del - tretji_del; 
                    end
                end
            end
    end
end

alfe =  A\b;
P_res = [P,zeros(n,1)];

for i = 1:m
    P_res(neRobni(i),3) = alfe(i);
end

for i = 1:length(Robni)
    P_res(Robni(i),3) = g(P_res(Robni(i),1), P_res(Robni(i),2));
end

u = triangulation(t.ConnectivityList,P_res);

end