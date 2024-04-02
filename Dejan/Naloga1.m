f1 = @(x) 0*x + 0.5;
f2 = @(x) 2*x;
f3 = @(x) 2 - 2*x;
f4 = @(x) 2*x - 1;
f5 = @(x) -2*x + 1;

X = linspace(0,1,101);
plot(X, f1(X), X, f2(X), X, f3(X), X, f4(X), X, f5(X)), ylim([0, 1])
pause
T = [1,4,5;1,2,5;2,5,6;2,3,6;3,6,7;4,5,8;5,8,9;5,6,9;6,9,10;6,7,10];    % Ogljisca trikotnikov
X = [0; 0.5; 1; 0; 0.25; 0.75; 1; 0; 0.5; 1];
Y = [0; 0; 0; 0.5; 0.5; 0.5; 0.5; 1; 1; 1];

TR = triangulation(T,X,Y);            % Definiramo triangulacijo
TR.ConnectivityList;                % Seznam povezav
TR.Points;                          % Seznam tock
triplot(TR)
pause
i = 1;
Z = zeros(1, 10);
Z(i) = 1;
TO = triangulation(T,X(:),Y(:),Z(:));
trimesh(T,X,Y,Z);
colormap([0 0 1])
title(sprintf('h_{%g}', i));
xlabel('x');
ylabel('y');
zlabel(sprintf('h_{%g}(x,y)', i))
pause
grid on;
for i=2:10
    subplot(3,3,i-1)
    Z = zeros(1, 10);
    Z(i) = 1;
    TO = triangulation(T,X(:),Y(:),Z(:));
    trimesh(T,X,Y,Z);
    colormap([0 0 1])
    title(sprintf('h_{%g}', i));
    xlabel('x');
    ylabel('y');
    zlabel(sprintf('h_{%g}(x,y)', i))
end
