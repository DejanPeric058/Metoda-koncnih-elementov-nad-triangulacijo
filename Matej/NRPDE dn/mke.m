function u = mke(p,q,r,f,g,T)
PT = T.Points;
CT = T.ConnectivityList;
st_T = size(CT,1);
n = length(PT);
FB = freeBoundary(T);

notranje = [];
robne = [];

for i = 1:n
    if ismember(i,FB)
        notranje = [notranje,i];
    else
        robne = [robne,i];
    end
end
    
m = length(notranje);
A = zeros(m,m);
b = zeros(m,1);

for i = 1:st_T
    T = CT(i,:);
    trikotnik = [PT(T(1),:);PT(T(2),:);PT(T(3),:)];
    
    for j = 1:3
        [je_notranja1,id1] = ismember(T(j),notranje);
        G1 = zeros(3,1);
        if je_notranja1
            G1(j) = 1;
        else
            if isequal(class(g),'sym')
                syms x y    
                G1(j) = subs(g,[x y],PT(T(j),:));
                clear x y
            else
                G1(j) = g(PT(T(j),1),PT(T(j),2));
            end
        end
        
        for k = i:3
            [je_notranja2,id2] = ismember(T(k),notranje);
            G2 = zeros(3,1);
            if je_notranja2
                    G2(k) = 1;
            else
                if isequal(class(g),'sym')
                    syms x y
                    G2(k) = subs(g,[x y],PT(T(k),:));
                    clear x y
                else
                    G2(k) = g(PT(T(k),1),PT(T(k),2));
                end
            end
            
            [phi1,phi_x1,phi_y1] = linFunc(trikotnik,G1);
            [phi2,phi_x2,phi_y2] = linFunc(trikotnik,G2);
            
            int_p = priblInteg(trikotnik(:,1),trikotnik(:,2),@(x,y) p(x,y));
            int_q = priblInteg(trikotnik(:,1),trikotnik(:,2),@(x,y) q(x,y));
            
            clen1 = int_p * phi_x1 * phi_x2;
            clen2 = int_q * phi_y1 * phi_y2;
            clen3 = priblInteg(trikotnik(:,1),trikotnik(:,2),@(x,y) r(x,y).*phi1(x,y).*phi2(x,y));
            
            if je_notranja1&&je_notranja2
                A(id1,id2) = A(id1,id2) + clen1 + clen2 + clen3;
                A(id2,id1) = A(id1,id2); 
                
                if j == k
                    [phi_0,~,~] = linFunc(trikotnik,G1);
                    clen4 = priblInteg(trikotnik(:,1),trikotnik(:,2),@(x,y) f(x,y).*phi_0(x,y));
                
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
 
alpha_u = A\b;
PT_u = [PT, zeros(n,1)];

for i = 1:m
    PT_u(notranje(i),3) = alpha_u(i);
end

for i = 1:(n-m)
    if isequal(class(g),'sym')
        syms x y
        PT_u(robne(i),3) = subs(g,[x y],PT_u(robne(i),:));
        clear x y
    else
        PT_u(robne(i),3) = g(PT_u(robne(i),1),PT_u(robne(i),2));
    end
end

u = triangulation(CT,PT_u);
end



