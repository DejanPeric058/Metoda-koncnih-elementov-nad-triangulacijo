function [phi,phi_x,phi_y] = linFunc(E,g)
A = [ones(3,1),E];
B = A\g;

phi = @(x,y) B(1) + B(2)*x + B(3)*y;
phi_x = B(1);
phi_y = B(2);
end

