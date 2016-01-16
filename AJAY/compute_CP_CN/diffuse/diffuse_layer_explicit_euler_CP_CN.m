
function [ Cpn Cnn ] = diffuse_layer_explicit_euler_CP_CN( P, S, D, F, t, Cpo, Cpn, Cno, Cnn, Un )

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
    
   [ Cpo Cno ] = diffuse_boundary_condition_CP_CN( P, S, D, F, t, Cpo, Cno, Un );
  
    for j = 2:Ny+1
    
        for i = 2:Nx+1
      
            if ( F.LAYER(j,i) == F.DIFFUSE )

                % CATION
                % CALCULATING TERM: P1
                P1 = dt*D_p*(  ((Cpo(j,i+1) - 2*Cpo(j,i) + Cpo(j,i-1))/(hx*hx)) ...
                             + ((Cpo(j+1,i) - 2*Cpo(j,i) + Cpo(j-1,i))/(hy*hy)) );
                         
                % CALCULATING TERM: P2 
%                 P2 = dt*m_p*(  ((Cpo(j,i+1) - Cpo(j,i-1))*(Un(j,i+1) - Un(j,i-1))/(4*hx*hx)) ...
%                              + ((Cpo(j+1,i) - Cpo(j-1,i))*(Un(j+1,i) - Un(j-1,i))/(4*hy*hy)) );
                P2 = dt * m_p *((Cpo(j,i) - Cpo(j-1,i))*(Un(j,i) - Un(j-1,i))/hx^2 + (Cpo(j,i) - Cpo(j,i-1))*(Un(j,i) - Un(j,i-1))/hy^2 );
                
                % CALCULATING TERM: P3
                P3 = dt*m_p*Cpo(j,i)*(  ((Un(j,i+1) - 2*Un(j,i) + Un(j,i-1))/(hx*hx)) ...
                                      + ((Un(j+1,i) - 2*Un(j,i) + Un(j-1,i))/(hy*hy)) );
        
                Cpn(j,i) = Cpo(j,i) + P1 + P2 + P3;
        
        
                % ANION
                % CALCULATING TERM: N1
                N1 = dt*D_n*(  ((Cno(j,i+1) - 2*Cno(j,i) + Cno(j,i-1))/(hx*hx)) ...
                             + ((Cno(j+1,i) - 2*Cno(j,i) + Cno(j-1,i))/(hy*hy)) );
                   
                % CALCULATING TERM: N2
%                 N2 = dt*m_n*(  ((Cno(j,i+1) - Cno(j,i-1))*(Un(j,i+1) - Un(j,i-1))/(4*hx*hx)) ...
%                              + ((Cno(j+1,i) - Cno(j-1,i))*(Un(j+1,i) - Un(j-1,i))/(4*hy*hy)) );
                    
                % CALCULATING TERM: N3
                N3 = dt*m_n*Cno(j,i)*(  ((Un(j,i+1) - 2*Un(j,i) + Un(j,i-1))/(hx*hx)) ...
                                      + ((Un(j+1,i) - 2*Un(j,i) + Un(j-1,i))/(hy*hy)) );
        
                Cnn(j,i) = Cno(j,i) + N1 - N3;         
                                                            
            end
    
        end

    end
                                   
end