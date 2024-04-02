%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Skiciranje linearnih funkcij
f1 = @(x) 0*x + 0.5;
f2 = @(x) 2*x;
f3 = @(x) 2 - 2*x;
f4 = @(x) 2*x - 1;
f5 = @(x) -2*x + 1;

X = linspace(0,1,101);
plot(X, f1(X), X, f2(X), X, f3(X), X, f4(X), X, f5(X)), ylim([0, 1])
pause
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Definiranje triangulacije
T1 = [1,4,5;1,2,5;2,5,6;2,3,6;3,6,7;4,5,8;5,8,9;5,6,9;6,9,10;6,7,10];    % Ogljisca trikotnikov
X1 = [0; 0.5; 1; 0; 0.25; 0.75; 1; 0; 0.5; 1];
Y1 = [0; 0; 0; 0.5; 0.5; 0.5; 0.5; 1; 1; 1];

TR1 = triangulation(T1,X1,Y1);              % Definiramo triangulacijo
triangles = TR1.ConnectivityList            % Seznam povezav
points = TR1.Points                         % Seznam tock
hold on
triplot(TR1)                                % Narišemo triangulacijo
hold off
pause
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Skiciranje h_1
i = 1;
Z1 = zeros(1, 10);
Z1(i) = 1;                                  % Vektor e1

trimesh(T1,X1,Y1,Z1);                       
colormap([0 0 1])
title(sprintf('h_{%g}', i));
xlabel('x');
ylabel('y');
zlabel(sprintf('h_{%g}(x,y)', i))
hold off
pause
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Skiciranje ostalih h_i

for i=2:10
    subplot(3,3,i-1)
    Z1 = zeros(1, 10);
    Z1(i) = 1;
    
    trimesh(T1,X1,Y1,Z1);
    colormap([0 0 1])
    title(sprintf('h_{%g}', i));
    xlabel('x');
    ylabel('y');
    zlabel(sprintf('h_{%g}(x,y)', i))
end
hold off
pause
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Prikaz triangulacije z vgajenim ukazom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% delaunay
subplot(1,1,1)
hold on
[X2, Y2] = meshgrid(linspace(0, 1, 6));     % Mreža 6x6 točk
X2 = X2(:);                                 % Iz matrike v vektor
Y2 = Y2(:);                                 % Iz matrike v vektor
T2 = delaunay(X2, Y2);                      % Array točk, ki so povezane
TR2 = triangulation(T2,X2,Y2);              % Definiramo triangulacijo
triplot(TR2)                                % Jo narišemo
hold off
pause
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Prikaz h_i na tej triangulaciji
Z2 = zeros(1,36);
Z2(16) = 1;                                 % Naredimo nek enotski vektor
trimesh(T2,X2,Y2,Z2(:));                    % Narišemo h_i
colormap([0 0 1])
title(sprintf('h_{%g}', i));
xlabel('x');
ylabel('y');
zlabel(sprintf('h_{%g}(x,y)', i))
hold off
