
x = (0:1/(P.Nx-1):1); y = (0:1/(P.Ny-1):1);
[X Y] = meshgrid(x,y);

figure(2)
hold on
clf

subplot(2,2,1)
surf(X,Y,Cpt([2:P.Nx+1],[2:P.Ny+1],i)/(P.NA*(P.Debye)^3))
% shading faceted; light; lighting flat;
shading interp; light; lighting phong;
xlabel('[Scaled L]')
ylabel('[Scaled L]')
zlabel('[mol m^{-3}]')

subplot(2,2,2)
surface(X,Y,Cpt([2:P.Nx+1],[2:P.Ny+1],i)/(P.NA*(P.Debye)^3))
% shading faceted; light; lighting flat;
shading interp; light; lighting phong;
xlabel('[Scaled L]')
ylabel('[Scaled L]')

subplot(2,2,[3,4])
contour(X,Y,Cpt([2:P.Nx+1],[2:P.Ny+1],i)/(P.NA*(P.Debye)^3))
xlabel('[Scaled L]')
ylabel('[Scaled L]')

title('CATION IN PNP MODEL')
colorbar

figure(3)
hold on
clf

subplot(2,2,1)
surf(X,Y,Cnt([2:P.Nx+1],[2:P.Ny+1],i)/(P.NA*(P.Debye)^3))
% shading faceted; light; lighting flat;
% shading interp; light; lighting phong;
xlabel('[Scaled L]')
ylabel('[Scaled L]')
zlabel('[mol m^{-3}]')

subplot(2,2,2)
surface(X,Y,Cnt([2:P.Nx+1],[2:P.Ny+1],i)/(P.NA*(P.Debye)^3))
% shading faceted; light; lighting flat;
% shading interp; light; lighting phong;
xlabel('[Scaled L]')
ylabel('[Scaled L]')

subplot(2,2,[3,4])
contour(X,Y,Cnt([2:P.Nx+1],[2:P.Ny+1],i)/(P.NA*(P.Debye)^3))
xlabel('[Scaled L]')
ylabel('[Scaled L]')

title('ANION IN PNP MODEL')
colorbar
