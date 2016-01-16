
function [ residual ] = diffuse_layer_sor_residual( P, S, D, F, t, Cpo, Cno, Uo, Un )

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
        
                % CALCULATING RHS
                RHS = c*(Cno(j,i) - Cpo(j,i));
        
                residual = residual + ( RHS  ...
                                         - (Un(j,i+1) - 2*Un(j,i) + Un(j,i-1))/(hx*hx) ...
                                         - (Un(j+1,i) - 2*Un(j,i) + Un(j-1,i))/(hy*hy) )^2;                               
                                      
            end
    
        end
    
    end
    
    residual = sqrt(residual/(Nx*Ny));
  
end