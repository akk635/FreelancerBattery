
tic;

[ P S D F Uo Un Ut Cpo Cpn Cpt Cno Cnn Cnt ] = init( filename );

% IMPOSE BOUNDARY CONDITIONS
t = 0;
[ Uo ]      = diffuse_boundary_condition_U    ( P, S, D, F, t, Uo );
[ Cpo Cno ] = diffuse_boundary_condition_CP_CN( P, S, D, F, t, Cpo, Cno, Uo );
iteration_t = 0; counter = 1;
% STORE RESULT AT t = 0
Ut(:,:,counter) = Uo; Cpt(:,:,counter) = Cpo; Cnt(:,:,counter) = Cno;

for t = D.t_start:D.dt:D.t_end
     
    [ Un Cpn Cnn iteration_DU residual_DU ] = diffuse_layer( P, S, D, F, t, Cpo, Cpn, Cno, Cnn, Uo, Un);
    
    [Csn] = particle_Diffsion(P, F, t, Csn, Cso, Cpn, i_app );
    
                    
    disp([                                             ...
          'Timestep : ', sprintf(num2str(iteration_t)),  ...
                         sprintf('\t \t'),               ...
          'Iter DU : ',  sprintf(num2str(iteration_DU)), ...
                         sprintf('\t \t'),               ...
          'Res DU : ',   sprintf(num2str(residual_DU)),  ...
                         sprintf('\t \t')
         ]);    
            
     % STORE RESULT AT t
     if( mod(iteration_t,5) == 0 ) 
        Ut(:,:,counter) = Un; Cpt(:,:,counter) = Cpn; Cnt(:,:,counter) = Cnn;
                counter = counter + 1;
     end
     
     Uo  = Un; Cpo = Cpn; Cno = Cnn;
     iteration_t = iteration_t + 1;        
     
end % END OF FOR
 
Total_Time = toc;
