
function [ residual ] = diffuse_layer_gs_residual_formulation_2( P, S, D, F, t, Cpo, Cno, Uo, Un )

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

    residual = 0;
  
    for j = 2:Ny+1
    
        for i = 2:Nx+1
      
            if ( F.LAYER(j,i) == F.DIFFUSE )
        
                % CALCULATING TERM: K1
                K1 = 1 + c*dt*(m_n*Cno(j,i) + m_p*Cpo(j,i));
  
                % CALCULATING TERM: K2
                K2 = c*dt*m_n*(  (Cno(j,i+1) - Cno(j,i-1))*(Uo(j,i+1) - Uo(j,i-1))/(4*hx*hx) ...
                               + (Cno(j+1,i) - Cno(j-1,i))*(Uo(j+1,i) - Uo(j-1,i))/(4*hy*hy) );
  
                % CALCULATING TERM: K3
                K3 = c*dt*m_p*(  (Cpo(j,i+1) - Cpo(j,i-1))*(Uo(j,i+1) - Uo(j,i-1))/(4*hx*hx) ...
                               + (Cpo(j+1,i) - Cpo(j-1,i))*(Uo(j+1,i) - Uo(j-1,i))/(4*hy*hy) );
  
                % CALCULATING TERM: K4
                K4 = c*(Cno(j,i) - Cpo(j,i)) + c*dt*(  D_n*(  (Cno(j,i+1) - 2*Cno(j,i) + Cno(j,i-1))/(hx*hx)   ...
                                                            + (Cno(j+1,i) - 2*Cno(j,i) + Cno(j-1,i))/(hy*hy) ) ...
                                                     - D_p*(  (Cpo(j,i+1) - 2*Cpo(j,i) + Cpo(j,i-1))/(hx*hx)   ...
                                                            + (Cpo(j+1,i) - 2*Cpo(j,i) + Cpo(j-1,i))/(hy*hy)) );
                
                RHS = (K4 - K2 - K3)/K1;
        
                residual = residual + ( RHS  ...
                                         - (Un(j,i+1) - 2*Un(j,i) + Un(j,i-1))/(hx*hx) ...
                                         - (Un(j+1,i) - 2*Un(j,i) + Un(j-1,i))/(hy*hy) )^2;                               
                                      
            end
    
        end
    
    end
    
    residual = sqrt(residual/(Nx*Ny));
  
end
