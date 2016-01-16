
x = (0:1/(P.Nx-1):1); 

figure(1)
hold on
clf
i=1;
subplot(2,1,1)
U2D = Ut([2:P.Ny+1],[2:P.Nx+1],i);
plot(x*P.Debye,U2D(8,:)*P.Ut)
xlabel('[Scaled L]')
ylabel('[V]')

subplot(2,1,2)
hold on
Cp2D = Cpt([2:P.Ny+1],[2:P.Nx+1],i);
Cn2D = Cnt([2:P.Ny+1],[2:P.Nx+1],i);
plot(x*P.Debye,Cp2D(8,:)/(P.NA*(P.Debye)^3))
plot(x*P.Debye,Cn2D(8,:)/(P.NA*(P.Debye)^3))
xlabel('[Scaled L]')
ylabel('[mol m^{-3}]')