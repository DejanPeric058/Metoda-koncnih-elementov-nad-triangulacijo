function a = funk_g(x,y)
if(0<=x)&&(x<=1)&&(y==0)
    a = x^3;
elseif(0<=x)&&(x<=1)&&(y==1)
    a = x^2;
elseif(x==0)&&(0<=y)&&(y<=1)
    a = sin(2*pi*y);
elseif(x==0)&&(0<=y)&&(y<=1)
    a = cos(2*pi*y);
else
    NaN
end