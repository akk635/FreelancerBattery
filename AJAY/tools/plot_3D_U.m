
x = (0:1/(P.Nx-1):1); y = (0:1/(P.Ny-1):1);
[X Y] = meshgrid(x,y);

figure(1)
hold on
clf

subplot(2,2,1)
surf(X*P.Debye,Y*P.Debye,Ut([2:P.Nx+1],[2:P.Ny+1],i)*P.Ut)
% shading faceted; light; lighting flat;
shading interp; light; lighting phong;
xlabel('[Scaled L]')
ylabel('[Scaled L]')
zlabel('[V]')

subplot(2,2,2)
surface(X,Y,Ut([2:P.Nx+1],[2:P.Ny+1],i)*P.Ut)
% shading faceted; light; lighting flat;
shading interp; light; lighting phong;
xlabel('[Scaled L]')
ylabel('[Scaled L]')

subplot(2,2,[3,4])
contour(X,Y,Ut([2:P.Nx+1],[2:P.Ny+1],i)*P.Ut)
xlabel('[Scaled L]')
ylabel('[Scaled L]')

title('SOLUTION FOR POTENTIAL IN PNP MODEL')
colorbar
