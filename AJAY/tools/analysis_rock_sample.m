
S = load(filename);

[Nx,Ny] = size(S.A);

% SPACE DISCRETIZATION
% GRID POINTS FOR X AXIS
hx = 1/(Nx-1);
% GRID POINTS FOR Y AXIS  
hy = 1/(Ny-1);
 
figure(10)
hold on;
title('GEOMETRICAL STRUCTURE OF ROCK SAMPLE 2D','FontSize',10);
x = 0:hx:1; y = 0:hy:1;
[X Y] = meshgrid(x,y);
surface(X,Y,S.A);

 