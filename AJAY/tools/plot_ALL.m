
  figure(1)
  potential = Ut([2:P.Nx+1],[2:P.Ny+1],:)*P.Ut;
  set(gca,'clim',[min(min(min(potential(potential>0)))) max(max(max(potential)))])
  surface(X,Y,Ut([2:P.Nx+1],[2:P.Ny+1],i)*P.Ut)
  shading interp;
  title('POTENTIAL')
  colorbar

  figure(2)
  cation = Cpt([2:P.Nx+1],[2:P.Ny+1],:)/(P.NA*(D.Debye)^3);
  set(gca,'clim',[min(min(min(cation(cation>0)))) max(max(max(cation)))])
  surface(X,Y,Cpt([2:P.Nx+1],[2:P.Ny+1],i)/(P.NA*(D.Debye)^3))
  shading interp; 
  title('CATION CONC.')
  colorbar
  
  figure(3)
  anion = Cnt([2:P.Nx+1],[2:P.Ny+1],:)/(P.NA*(D.Debye)^3);
  set(gca,'clim',[min(min(min(anion(anion>0)))) max(max(max(anion)))])
  surface(X,Y,Cnt([2:P.Nx+1],[2:P.Ny+1],i)/(P.NA*(D.Debye)^3))
  shading interp;
  title('ANION CONC.') 
  colorbar