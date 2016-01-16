
function [ Un ] = diffuse_layer_gs_solver_formulation_1( P, S, D, F, t, Cpo, Cno, Uo, Un )

    c     = D.c;
    sigma = D.sigma;
    D_p   = D.D_p;
    D_n   = D.D_n;
    m_p   = D.m_p;
    m_n   = D.m_n;
  
    Nx  = P.Nx;
    Ny  = P.Ny;
    hx  = D.hx;
    hy  = D.hy;
    hxy = (2/(hx*hx)) + (2/(hy*hy));  
    dt  = D.dt;
    
    Un = Uo;
  
    for j = 2:Ny+1
    
        for i = 2:Nx+1
      
            if ( F.LAYER(j,i) == F.DIFFUSE )
        
                % CALCULATING TERM: K1
                K1 = hxy;
                
                % CALCULATING TERM: K2
                K2 = c*(Cno(j,i) - Cpo(j,i));
  
                % CALCULATING TERM: K3
                K3 = (Un(j,i+1) + Un(j,i-1))/(hx*hx) +  (Un(j+1,i) + Un(j-1,i))/(hy*hy);
                
                Un(j,i) = (K3 - K2)/K1;
      
            end
    
        end

    end
    
    % IMPOSE BOUNDARY CONDITION
    Un = diffuse_boundary_condition_U( P, S, D, F, t, Un ); 
                                   
end
