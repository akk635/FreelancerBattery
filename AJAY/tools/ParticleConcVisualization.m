% visualize particle concentration in transient time steps
x = (0:1/(P.Nx-1):1); y = (0:1/(P.Ny-1):1);
[X Y] = meshgrid(x,y);

filename1 = 'potentialnConcentration.gif';
filename2 = 'solidConcentration.gif';
filename3 = 'negativeConcentration.gif';
for i = 1 : 100 : counter
      figure (2);
      curr_Csn = Cst([2:end-1],[2:end-1],i);
      set(gca,'clim',[P.Cs0 max(max(max(curr_Csn)))]);
      surface(X,Y,curr_Csn);
%       hold on;
%       contour(X,Y, curr_Csn);
      colorbar

%       figure (1);
%       subplot(1,2,1)
%       potential = Ut([2:P.Nx+1],[2:P.Ny+1],i);
%       set(gca,'clim',[min(min(min((potential)))) max(max(max(potential)))])
%       surface(X,Y,potential);
%       hold on;
%       contour(X,Y,potential);
%       title('Potential')
%       colorbar
%         
%       cation = Cpt([2:P.Nx+1],[2:P.Nx+1],i);
%       anion  = Cnt([2:P.Nx+1],[2:P.Ny+1],i);
%       MIN = [min(min(min(cation(cation>0)))), min(min(min(anion(anion>0))))];
%       MAX = [max(max(max(cation))), max(max(max(anion)))];
% 
%         subplot(1,2,2)
%         set(gca,'clim',[MIN(1) MAX(1)])
%           surface(X,Y,cation)
%           hold on;
%           contour(X,Y,cation);
%           title('Cation conc.')
%           colorbar
% 
%           figure(3);
%           set(gca,'clim',[MIN(2) MAX(2)])
%           surface(X,Y, anion)
%           hold on;
%           contour(X,Y,anion);
%           title('Anion conc.')
%           colorbar
          
                drawnow
      frame = getframe(2);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if i == 1;
          imwrite(imind,cm,filename2,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename2,'gif','WriteMode','append');
      end
%       frame = getframe(1);
%       im = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);
%       if i == 1;
%           imwrite(imind,cm,filename1,'gif', 'Loopcount',inf);
%       else
%           imwrite(imind,cm,filename1,'gif','WriteMode','append');
%       end
%             frame = getframe(3);
%       im = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);
%       if i == 1;
%           imwrite(imind,cm,filename3,'gif', 'Loopcount',inf);
%       else
%           imwrite(imind,cm,filename3,'gif','WriteMode','append');
%       end
          
          

          

end