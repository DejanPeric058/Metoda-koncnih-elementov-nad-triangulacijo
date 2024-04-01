function [phi,phi_x,phi_y] = linFunc(d,v)
% 

A = [v(1) v(2) 1; d(1,:) 1; d(2,:) 1];
D = A(2,1)*A(3,2)-A(3,1)*A(2,2)-A(1,1)*A(3,2)+A(3,1)*A(1,2)+A(1,1)*A(2,2)-A(2,1)*A(1,2);
a = (A(2,2)-A(3,2))/D; b = (A(3,1)-A(2,1))/D; c = (A(2,1)*A(3,2)-A(3,1)*A(2,3))/D;

phi = @(x,y) a*x + b*y + c;
phi_x = a;
phi_y = b;
end

