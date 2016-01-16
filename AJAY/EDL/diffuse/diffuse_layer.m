
function [ Un Cpn Cnn iteration_U residual_U ] = diffuse_layer( P, S, D, F, t, Cpo, Cpn, Cno, Cnn, Uo, Un )
      
    switch ( P.ALGO )
       
        case ( 1 ) % ALGO 1: Simple  GS Solver (U) + Explicit Euler (Cp Cn)    
            [ Un iteration_U residual_U ] = diffuse_layer_compute_U( P, S, D, F, t, Cpo, Cno, Uo, Un );
            [ Cpn Cnn ] = diffuse_layer_explicit_euler_CP_CN( P, S, D, F, t, Cpo, Cpn, Cno, Cnn, Un );   

       
        case ( 2 ) % ALGO 2: Simple SOR Solver (U) + Explicit Euler (Cp Cn)        
            [ Un iteration_U residual_U ] = diffuse_layer_compute_U( P, S, D, F, t, Cpo, Cno, Uo, Un );
            [ Cpn Cnn ] = diffuse_layer_explicit_euler_CP_CN( P, S, D, F, t, Cpo, Cpn, Cno, Cnn, Un );   

       
        case ( 3 ) % ALGO 3: Complex GS Solver (U) + Explicit Euler (Cp Cn)        
            [ Un iteration_U residual_U ] = diffuse_layer_compute_U( P, S, D, F, t, Cpo, Cno, Uo, Un );
            [ Cpn Cnn ] = diffuse_layer_explicit_euler_CP_CN( P, S, D, F, t, Cpo, Cpn, Cno, Cnn, Un );   
                
    end
    
end
    