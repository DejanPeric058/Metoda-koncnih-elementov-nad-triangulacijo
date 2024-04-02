function integral = triintegral(k,A)
phi1 = @(u,v) u*A(1,1)+v*A(2,1)+(1-u-v)*A(3,1);
phi2 = @(u,v) u*A(1,2)+v*A(2,2)+(1-u-v)*A(3,2);
det = (A(1,1) - A(3,1))*(A(2,2) - A(3,2)) - (A(2,1) - A(3,1))*(A(1,2) - A(3,2));
integral = integral2(@(u,v)k(phi1(u,v),phi2(u,v))*abs(det),0,1,0,@(u)1-u);
end

