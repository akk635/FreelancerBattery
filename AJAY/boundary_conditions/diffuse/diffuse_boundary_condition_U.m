
function [ Uo ] = diffuse_boundary_condition_U( P, S, D, F, t, Uo )
    
    c     = D.c;
    sigma = D.sigma;
    D_p   = D.D_p;
    D_n   = D.D_n;
    m_p   = D.m_p;
    m_n   = D.m_n;
  
    Nx = P.Nx;
    Ny = P.Ny;
    hx = D.hx;
    hy = D.hy;
    dt = D.dt;
  
    %
    %   A ----------------------- D
    %   |                         |
    %   |                         |
    %   |                         |
    %   B ----------------------- C
    %    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXTERNAL BOUNDARY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (i.e. SIDE AB - APPLIED AC OR DC CURRENT)
    
    % WEST: DIRICHLET 
    Uo(:,1) = P.U_ext/P.Ut;         
    
    % EAST: NEUMANN + dirichlet zero potential (Robin boundary cdn)
    Uo(:,Nx+2) = 0;
    
    % NORTH: NEUMANN
    Uo(1,:) = Uo(2,:);

    % SOUTH: NEUMANN
    Uo(Ny+2,:) = Uo(Ny+1,:); 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INTERNAL BOUNDARY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (i.e. ROCK-FLUID INTERFACE) % considering electrostatic boundary
    % condition
    
    for j = 2:Ny+1
    
        for i = 2:Nx+1 
                        
            % DIRICHLET: INTERNAL BOUNDARY
            if ( (F.LAYER(j,i) == F.OBJECT) & (P.BC_U == 1) )
                
                Uo(j,i) = sigma;
              
            end            
            
            
            % NEUMANN ZERO: INTERNAL BOUNDARY
            if ( (F.LAYER(j,i) == F.OBJECT) & (P.BC_U == 2) )
                
%                 Uo(j,i) = sigma;                   
                                 
                switch ( F.FLAG(j,i) )
                         
                    % 1 DIRECTIONS: TOTAL 4 COMBINATIONS                   

                    case ( F.B_E )  % EAST BOUNDARY
                        
                        Uo(j,i) = Uo(j,i+1);
                                              
                    case ( F.B_W )  % WEST BOUNDARY
                        
                        Uo(j,i) = Uo(j,i-1);                                         
                                               
                    case ( F.B_N )  % NORTH BOUNDARY
                        
                        Uo(j,i) = Uo(j-1,i);                    
                
                    case ( F.B_S )  % SOUTH BOUNDARY
                         
                        Uo(j,i) = Uo(j+1,i); 
                        
                    % 2 DIRECTIONS: TOTAL 6 COMBINATIONS
                    
                    case ( F.B_N + F.B_E )  % NORTH EAST BOUNDARY
 
                        Uo(j,i) = (  Uo(j-1,i)...
                                   + Uo(j,i+1) )/2; 
                            
                    case ( F.B_N + F.B_W )  % NORTH WEST BOUNDARY
 
                        Uo(j,i) = (  Uo(j-1,i) ...
                                   + Uo(j,i-1) )/2;                     
                     
                    case ( F.B_S + F.B_E )  % SOUTH EAST BOUNDARY

                        Uo(j,i) = (  Uo(j+1,i) ...
                                   + Uo(j,i+1) )/2; 
                                                                                               
                    case ( F.B_S + F.B_W )  % SOUTH WEST BOUNDARY

                        Uo(j,i) = (  Uo(j+1,i) ...
                                   + Uo(j,i-1) )/2; 
                                       
                    case ( F.B_N + F.B_S )  % NORTH SOUTH BOUNDARY    
                        
                        Uo(j,i) = (  Uo(j-1,i) ...                       
                                   + Uo(j+1,i) )/2;                       
                                                                      
                    case ( F.B_E + F.B_W )  % EAST WEST BOUNDARY    
                        
                        Uo(j,i) = (  Uo(j,i-1) ...                         
                                  +  Uo(j,i+1) )/2;                         
                                               
                    % 3 DIRECTIONS: TOTAL 4 COMBINATIONS   
                                                
                    case ( F.B_N + F.B_E + F.B_W )  % NORTH EAST WEST BOUNDARY

                        Uo(j,i) = (  Uo(j-1,i) ...
                                   + Uo(j,i+1) ...
                                   + Uo(j,i-1) )/3; 
     
                    case ( F.B_S + F.B_E + F.B_W )  % SOUTH EAST WEST BOUNDARY

                        Uo(j,i) = (  Uo(j+1,i) ...
                                   + Uo(j,i+1) ...
                                   + Uo(j,i-1) )/3;                        
                                                                           
                    case ( F.B_N + F.B_S + F.B_W )  % NORTH SOUTH WEST BOUNDARY

                        Uo(j,i) = (  Uo(j-1,i) ...
                                   + Uo(j+1,i) ...
                                   + Uo(j,i-1) )/3; 
                        
                    case ( F.B_N + F.B_S + F.B_E )  % NORTH SOUTH EAST BOUNDARY

                        Uo(j,i) = (  Uo(j-1,i) ...
                                   + Uo(j+1,i) ...
                                   + Uo(j,i+1) )/3; 
                                                
                    otherwise
                        
                end % END OF SWITCH
              
            end

            
            % NEUMANN UNZERO: INTERNAL BOUNDARY
            if ( (F.LAYER(j,i) == F.OBJECT) & (P.BC_U == 3) )
                
               % ASSIGN A SURFACE CHARGE DENSITY
                surface_charge_density = P.sigma * 1e-3; 
                RHS = (surface_charge_density/(P.epsilonO))*(P.Debye/P.Ut);
                                                                
                switch ( F.FLAG(j,i) )
                         
                    % 1 DIRECTIONS: TOTAL 4 COMBINATIONS                   

                    case ( F.B_E )  % EAST BOUNDARY
                        
                        Uo(j,i) = Uo(j,i+1) + RHS*hx;
                                              
                    case ( F.B_W )  % WEST BOUNDARY
                        
                        Uo(j,i) = Uo(j,i-1) + RHS*hx;                                        
                                               
                    case ( F.B_N )  % NORTH BOUNDARY
                        
                        Uo(j,i) = Uo(j-1,i) + RHS*hy;                  
                
                    case ( F.B_S )  % SOUTH BOUNDARY
                         
                        Uo(j,i) = Uo(j+1,i) + RHS*hy;
                        
                    % 2 DIRECTIONS: TOTAL 6 COMBINATIONS
                    
                    case ( F.B_N + F.B_E )  % NORTH EAST BOUNDARY
 
                        Uo(j,i) = (  Uo(j-1,i) + RHS*hx ...
                                   + Uo(j,i+1) + RHS*hy )/2; 
                            
                    case ( F.B_N + F.B_W )  % NORTH WEST BOUNDARY
 
                        Uo(j,i) = (  Uo(j-1,i) + RHS*hx ...
                                   + Uo(j,i-1) + RHS*hy )/2;                     
                     
                    case ( F.B_S + F.B_E )  % SOUTH EAST BOUNDARY

                        Uo(j,i) = (  Uo(j+1,i) + RHS*hx ...
                                   + Uo(j,i+1) + RHS*hy )/2; 
                                                                                               
                    case ( F.B_S + F.B_W )  % SOUTH WEST BOUNDARY

                        Uo(j,i) = (  Uo(j+1,i) + RHS*hx ...
                                   + Uo(j,i-1) + RHS*hy )/2; 
                                       
                    case ( F.B_N + F.B_S )  % NORTH SOUTH BOUNDARY    
                        
                        Uo(j,i) = (  Uo(j-1,i) + RHS*hx ...                       
                                   + Uo(j+1,i) + RHS*hy )/2;                       
                                                                      
                    case ( F.B_E + F.B_W )  % EAST WEST BOUNDARY    
                        
                        Uo(j,i) = (  Uo(j,i-1) + RHS*hx ...                         
                                  +  Uo(j,i+1) + RHS*hy )/2;                         
                                               
                    % 3 DIRECTIONS: TOTAL 4 COMBINATIONS   
                                                
                    case ( F.B_N + F.B_E + F.B_W )  % NORTH EAST WEST BOUNDARY

                        Uo(j,i) = (  Uo(j-1,i) + RHS*hy ...
                                   + Uo(j,i+1) + RHS*hx ...
                                   + Uo(j,i-1) + RHS*hx )/3; 
     
                    case ( F.B_S + F.B_E + F.B_W )  % SOUTH EAST WEST BOUNDARY

                        Uo(j,i) = (  Uo(j+1,i) + RHS*hy ...
                                   + Uo(j,i+1) + RHS*hx ...
                                   + Uo(j,i-1) + RHS*hx )/3;                        
                                                                           
                    case ( F.B_N + F.B_S + F.B_W )  % NORTH SOUTH WEST BOUNDARY

                        Uo(j,i) = (  Uo(j-1,i) + RHS*hy ...
                                   + Uo(j+1,i) + RHS*hy ...
                                   + Uo(j,i-1) + RHS*hx )/3; 
                        
                    case ( F.B_N + F.B_S + F.B_E )  % NORTH SOUTH EAST BOUNDARY

                        Uo(j,i) = (  Uo(j-1,i) + RHS*hy  ...
                                   + Uo(j+1,i) + RHS*hy  ...
                                   + Uo(j,i+1) + RHS*hx )/3; 
                                                
                    otherwise
                        
                end % END OF SWITCH
              
            end            
        end

    end
    
end
