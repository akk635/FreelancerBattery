
x = (0:1/(P.Nx+1):1); y = (0:1/(P.Ny+1):1);
[X Y] = meshgrid(x,y);

D = size(Ut);

message  = 'Enter the skip step? ';
skip     = input(message)

D(3) = counter;
 for(i=1:skip:D(3)-1)

  figure(4)
  hold on
  
  clf
  
  subplot(2,2,1)
  mesh(X,Y,Ut(:,:,i))
 
  subplot(2,2,2)
  surface(X,Y,Ut(:,:,i))
  
  subplot(2,2,[3,4])
  contour(X,Y,Ut(:,:,i))
  title('SOLUTION FOR POTENTIAL IN PNP MODEL') 
  colorbar
  
  pause(1)
 end 
