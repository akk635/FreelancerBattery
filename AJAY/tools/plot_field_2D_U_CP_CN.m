
x = (0:1/(P.Nx-1):1); y = (0:1/(P.Ny-1):1);
[X Y] = meshgrid(x,y);

hx = 1/(P.Nx-1);
hy = 1/(P.Ny-1);

D3 = size(Ut);

message  = 'Skip Time Point: Enter the skip step? ';
skip     = input(message);

for(i=1:skip:D3(3)-1)
    
  figure(6)
  hold on
  
  clf

  % SCALING CATION AND ANION CONCENTRATIONS
  % cation = Cpt([2:P.Nx+1],[2:P.Ny+1],end);
  % anion  = Cnt([2:P.Nx+1],[2:P.Ny+1],end);
  % MIN = min(min(min(min(cation(cation>0)))),min(min(min(anion(anion>0)))));
  % MAX = max(max(max(max(cation))),max(max(max(anion)))); 
  % set(gca,'clim',[MIN MAX])  
  
  % PLOT GRADIENT POTENTIAL
  subplot(2,2,2) 
  UV = Ut([2:P.Nx+1],[2:P.Ny+1],i);

  [ux,uy] = gradient(UV,hx,hy);

  contour(X,Y,UV)
  hold on
  quiver(x,y,ux,uy)
  colorbar
  
  % PLOT GRADIENT CATION CONCENTRATION
  subplot(2,2,3) 
  CpV = Cpt([2:P.Nx+1],[2:P.Ny+1],i);

  [cpx,cpy] = gradient(CpV,hx,hy);

  contour(X,Y,CpV)
  hold on
  quiver(x,y,cpx,cpy)
  colorbar  
  
  % PLOT GRADIENT ANION CONCENTRATION
  subplot(2,2,4) 
  CnV = Cnt([2:P.Nx+1],[2:P.Ny+1],i);

  [cnx,cny] = gradient(CnV,hx,hy);

  contour(X,Y,CnV)
  hold on
  quiver(x,y,cnx,cny)
  colorbar    
  
  % TIME STEP
  i
  pause(1)

end
