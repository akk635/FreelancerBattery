
function [ Un iteration_U residual_U ] = diffuse_layer_compute_U( P, S, D, F, t, Cpo, Cno, Uo, Un )

    Nx        = P.Nx;
    Ny        = P.Ny;
    max_iter  = P.max_iter;
    threshold = P.threshold;
  
    residual_U = 1; iteration_U = 0;
    
    switch ( P.ALGO )
        
        % ALGO 1: Simple  GS Solver (U) + Explicit Euler (Cp Cn)
        % ALGO 2: Simple SOR Solver (U) + Explicit Euler (Cp Cn)
        % ALGO 3: Complex GS Solver (U) + Explicit Euler (Cp Cn)    
        
        case ( 1 )
  
            while( (residual_U > threshold) & (iteration_U < max_iter) )
                Un          = diffuse_layer_gs_solver_formulation_1  ( P, S, D, F, t, Cpo, Cno, Uo, Un );
                residual_U  = diffuse_layer_gs_residual_formulation_1( P, S, D, F, t, Cpo, Cno, Uo, Un );
                iteration_U = iteration_U + 1;
                Uo = Un;    
            end % END OF WHILE
            Un = Uo;
            
        case ( 2 )
  
            while( (residual_U > threshold) & (iteration_U < max_iter) )
                Un          = diffuse_layer_sor_solver  ( P, S, D, F, t, Cpo, Cno, Uo, Un );
                residual_U  = diffuse_layer_sor_residual( P, S, D, F, t, Cpo, Cno, Uo, Un );
                iteration_U = iteration_U + 1;
                Uo = Un;    
            end % END OF WHILE
%             Un = Uo;
            
        case ( 3 )
  
            while( (residual_U > threshold) & (iteration_U < max_iter) )
                Un          = diffuse_layer_gs_solver_formulation_2  ( P, S, D, F, t, Cpo, Cno, Uo, Un );
                residual_U  = diffuse_layer_gs_residual_formulation_2( P, S, D, F, t, Cpo, Cno, Uo, Un );
                iteration_U = iteration_U + 1;
                Uo = Un;    
            end % END OF WHILE
            Un = Uo;            
            
    end % END OF SWITCH            
    
end
