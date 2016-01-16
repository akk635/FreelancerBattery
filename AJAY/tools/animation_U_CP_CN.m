
x = (0:1/(P.Nx-1):1); y = (0:1/(P.Ny-1):1);
[X Y] = meshgrid(x,y);

hx = 1/(P.Nx-1);
hy = 1/(P.Ny-1);

D3 = size(Ut);

message  = 'Skip Time Point: Enter the skip step? ';
skip     = input(message);

message = 'Choice of Visualization: 1 for Contour OR 2 for Surface? ';
visual  = input(message);

D3(3) = counter;

for (i=1:skip:D3(3)-1)
  
  figure(5)
  hold on
  
  clf

  % PLOT POTENTIAL
  subplot(2,2,2) 
  potential = Ut([2:P.Nx+1],[2:P.Ny+1],i);
  set(gca,'clim',[min(min(min((potential)))) max(max(max(potential)))])
  if (visual == 2)
    surface(X,Y,Ut([2:P.Nx+1],[2:P.Ny+1],i))
  else
    contour(X,Y,Ut([2:P.Nx+1],[2:P.Ny+1],i))
    %[ux,uy] = gradient(Ut([2:P.Nx+1],[2:P.Ny+1],i),hx,hy);
    %quiver(x,y,ux,uy)
  end
  %shading interp;
  title('POTENTIAL')
  colorbar

  % SCALING CATION AND ANION CONCENTRATIONS
  cation = Cpt(:,:,i);
  anion  = Cnt([2:P.Nx+1],[2:P.Ny+1],i);
  MIN = [min(min(min(cation(cation>0)))),min(min(min(anion(anion>0))))];
  MAX = [max(max(max(cation))),max(max(max(anion)))];
  
  % PLOT CATION CONCENTRATION
  subplot(2,2,3)

  set(gca,'clim',[MIN(1) MAX(1)])
  if(visual == 2) 
    surface(X,Y,Cpt([2:P.Nx+1],[2:P.Ny+1],i))
  else 
    contour(X,Y,Cpt([2:P.Nx+1],[2:P.Ny+1],i))  
    %[cpx,cpy] = gradient(Cpt([2:P.Nx+1],[2:P.Ny+1],i),hx,hy);
    %quiver(x,y,cpx,cpy)    
  end
  %shading interp; 
  title('CATION CONC.')
  colorbar
  
  % PLOT ANION CONCENTRAION
  subplot(2,2,4)
 
  set(gca,'clim',[MIN(2) MAX(2)])
  if(visual == 2)
    surface(X,Y,Cnt([2:P.Nx+1],[2:P.Ny+1],i))
  else
    contour(X,Y,Cnt([2:P.Nx+1],[2:P.Ny+1],i))  
    %[cnx,cny] = gradient(Cnt([2:P.Nx+1],[2:P.Ny+1],i),hx,hy);
    %quiver(x,y,cnx,cny)       
  end
  %shading interp;
  title('ANION CONC.') 
  colorbar
  
  % TIME STEP
  i
  pause(0.5)
  
end 
