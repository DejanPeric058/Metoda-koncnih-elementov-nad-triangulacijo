points = [0 0; 0 0.5; 0 1; 0.25 0.5; 0.5 0; 0.5 1; 0.75 0.5; 1 0; 1 0.5; 1 1];
connect = [1 2 4; 2 3 4; 1 4 5; 3 4 6; 4 5 7; 4 6 7; 5 7 8; 6 7 10; 7 8 9; 7 9 10];
T = triangulation(connect, points);
[A,B]=freeBoundary(T);

%hold on
%triplot(T)
%plot(B(:,1),B(:,2),'o')
%hold off

points2 = [points zeros(10,1); 0.25 0.5 1; 0.75 0.5 1];
connect2 = [1 2 11; 2 3 11; 1 4 11; 3 4 11; 4 5 11; 4 6 11; 4 5 12; 4 6 12; 5 7 12; 6 7 12; 7 8 12; 7 9 12;
1 4 11; 2 4 11; 1 5 11; 3 6 11; 4 7 11; 4 7 11; 4 7 12; 4 7 12; 5 8 12; 6 10 12; 7 9 12; 7 10 12;
2 4 11; 3 4 11; 4 5 11; 4 6 11; 5 7 11; 6 7 11; 5 7 12; 6 7 12; 7 8 12; 7 10 12; 8 9 12; 9 10 12];

T2 = triangulation(connect2, points2);
[A,B]=freeBoundary(T2);


trisurf(connect2, points2(:,1), points2(:,2), points2(:,3))


%[f,fx,fy] = linFunc([0 0; 1/2 0],[1/4 1/2]);
