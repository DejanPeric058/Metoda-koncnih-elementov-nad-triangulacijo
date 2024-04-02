function g = fun_g(x,y)
    if y == 0
        g = x^3;
    elseif y == 1
        g = x^2;
    elseif x == 0
        g = sin(2*pi*y);
    elseif x == 1
        g = cos(2*pi*y);
    else
        g = 0;
    end
end