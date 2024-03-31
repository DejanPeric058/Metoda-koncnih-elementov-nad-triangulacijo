function odgovor = fun_g(x,y)
if y == 0
    odgovor = x^3;
elseif y == 1
    odgovor = x^2;
elseif x == 0
    odgovor = sin(2*pi*y);
elseif x == 1
    odgovor = cos(2*pi*y);
else
    odgovor = 0;
end
end