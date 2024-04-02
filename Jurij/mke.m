function u = mke(tri, p, q, r, f, g)
% Opis:
%  mke izracuna priblizek za resitev parcialne diferencialne enacbe
%   - d/dx (p(x,y) du/dx) - d/dy (q(x,y) du/dy) + r(x,y) u = f(x,y)
%  z robnim pogojem u = g po metodi koncnih elementov z zveznimi odsekoma
%  linearnimi funkcijami nad triangulacijo
%
% Definicija
%  u = mke(tri,p,q,r,f,g)
%
% Vhodni podatki:
%  tri      triangulacija obmocja, predstavljena z razredom triangulation,
%  p,q,r,f  funkcije dveh spremenljivk, ki dolocajo parcialno diferencialno
%           enacbo,
%  g        funkcija dveh spremenljivk, ki doloca vrednost resitve na robu
%           obmocja
%
% Izhodni podatek:
%  u        3D triangulacija, predstavljena z razredom triangulation, ki
%           doloca zvezno odsekoma linearno funkcijo, ki je priblizek za
%           resitev robnega problema po metodi najmanjsih kvadratov
    num_triangles = size(tri.ConnectivityList, 1);
    
    all_vertices = tri.Points;
    all_vertices_indices = (1:size(tri.Points, 1))';
    
    boundary_indices = unique(freeBoundary(tri));
    boundary_vertices = all_vertices(boundary_indices, :);

    interior_indices = setdiff(all_vertices_indices, boundary_indices);
    interior_vertices = all_vertices(interior_indices, :);
    
    n = length(all_vertices);
    m = length(interior_vertices);
    A = zeros(m);
    b = zeros(m, 1);
    
    for k = 1:num_triangles
        current_triangle_indices = tri.ConnectivityList(k, :);
        current_triangle = tri.Points(current_triangle_indices, :);
        for i = 1:3
            [is_interior_i, ind_i] = ismember(current_triangle_indices(i), interior_indices);
            hi = zeros(3,1);
            if is_interior_i
                hi(i) = 1;
            else
                hi(i) = g(current_triangle(i, 1), current_triangle(i, 2));
            end
            for j = i:3
                [is_interior_j, ind_j] = ismember(current_triangle_indices(j), interior_indices);
                hj = zeros(3,1);
                if is_interior_j
                    hj(j) = 1;
                else
                    hj(j) = g(current_triangle(j, 1), current_triangle(j, 2));
                end
                
                [phi_i, phi_i_dx, phi_i_dy] = trilin(current_triangle, hi);
                [phi_j, phi_j_dx, phi_j_dy] = trilin(current_triangle, hj);
                
                p_integral = triintegral(@(x,y)p(x,y).*phi_i_dx.*phi_j_dx, current_triangle);
                q_integral = triintegral(@(x,y)q(x,y).*phi_i_dy.*phi_j_dy, current_triangle);
                r_integral = triintegral(@(x,y)r(x,y).*phi_j(x,y).*phi_i(x,y), current_triangle);
                total_integral = p_integral + q_integral + r_integral;

                if is_interior_i + is_interior_j == 2
                    A(ind_i,ind_j) = A(ind_i,ind_j) + total_integral;
                    A(ind_j,ind_i) = A(ind_i,ind_j);
                    if i == j
                        f_integral = triintegral(@(x,y) f(x,y).*phi_i(x,y), current_triangle);
                        b(ind_i) = b(ind_i) + f_integral;
                    end
                else
                    if is_interior_i
                        b(ind_i) = b(ind_i) - total_integral;
                    elseif is_interior_j
                        b(ind_j) = b(ind_j) - total_integral; 
                    end
                end
            end
        end
    end
    alfe =  A\b;
    P_tri = [tri.Points, zeros(n, 1)];

    for i = 1:length(interior_vertices)
        P_tri(interior_indices(i), 3) = alfe(i);
    end
    for i = 1:length(boundary_vertices)
        x = P_tri(boundary_indices(i), 1);
        y = P_tri(boundary_indices(i), 2);
        P_tri(boundary_indices(i), 3) = g(x, y);
    end
    u = triangulation(tri.ConnectivityList, P_tri);
end
