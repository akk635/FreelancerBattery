
function [ Cpo Cno ] = diffuse_boundary_condition_CP_CN( P, S, D, F, t, Cpo, Cno, Un )
    
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
    Uo = Un;
    
    if ( P.STERN == true )
        BOUNDARY = F.STERN;
    else
        BOUNDARY = F.OBJECT;
    end    
    
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
    
    % WEST electrostatic (no current at the bnd)
    Cpo(:,1) = (D_p + m_p*(Un(:,2) - Un(:,1))/2).*Cpo(:,2)./(D_p - m_p*(Un(:,2) - Un(:,1))/2);   
    Cno(:,1) = (D_n - m_n*(Un(:,2) - Un(:,1))/2).*Cno(:,2)./(D_n + m_n*(Un(:,2) - Un(:,1))/2);     

    % EAST electrostatic 
    Cpo(:,Nx+2) = (D_p - m_p*(Un(:,Nx+2) - Un(:,Nx+1))/2).*Cpo(:,Nx+1)./(D_p + m_p*(Un(:,Nx+2) - Un(:,Nx+1))/2);
    Cno(:,Nx+2) = (D_n + m_n*(Un(:,Nx+2) - Un(:,Nx+1))/2).*Cno(:,Nx+1)./(D_n - m_n*(Un(:,Nx+2) - Un(:,Nx+1))/2);
     
    % NORTH:
    Cpo(1,:) = Cpo(2,:);    
    Cno(1,:) = Cno(2,:);

    % SOUTH BOUNDARY
    Cpo(Ny+2,:) = Cpo(Ny+1,:);
    Cno(Ny+2,:) = Cno(Ny+1,:);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INTERNAL BOUNDARY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (i.e. ROCK-FLUID INTERFACE)
    
    for j = 2:Ny+1
        
        for i = 2:Nx+1 
            

            if ( (F.LAYER(j,i) == F.OBJECT) & (P.BC_CONC == 1) )
    
                switch ( F.FLAG(j,i) )
                         
                    % 1 DIRECTIONS: TOTAL 4 COMBINATIONS
                    
                    case ( F.B_E )  % EAST BOUNDARY
                         
                        Cpo(j,i) = Cpo(j,i+1) * ( 1 + ((m_p/D_p)*(Uo(j,i+1)-sigma)) );
                        
                        Cno(j,i) = Cno(j,i+1) * ( 1 - ((m_n/D_n)*(Uo(j,i+1)-sigma)) );
     
                    case ( F.B_W )  % WEST BOUNDARY
                         
                        Cpo(j,i) = Cpo(j,i-1) * ( 1 + ((m_p/D_p)*(Uo(j,i-1)-sigma)) );
                        
                        Cno(j,i) = Cno(j,i-1) * ( 1 - ((m_n/D_n)*(Uo(j,i-1)-sigma)) );
                       
                    case ( F.B_N )  % NORTH BOUNDARY
 
                        Cpo(j,i) = Cpo(j-1,i) * ( 1 + ((m_p/D_p)*(Uo(j-1,i)-sigma)) );
                        
                        Cno(j,i) = Cno(j-1,i) * ( 1 - ((m_n/D_n)*(Uo(j-1,i)-sigma)) );                         
               
                    case ( F.B_S )  % SOUTH BOUNDARY
                         
                        Cpo(j,i) = Cpo(j+1,i) * ( 1 + ((m_p/D_p)*(Uo(j+1,i)-sigma)) );
                        
                        Cno(j+1,i) = Cno(j+1,i) * ( 1 - ((m_n/D_n)*(Uo(j+1,i)-sigma)) ); 
                        
                        
                    % 2 DIRECTIONS: TOTAL 6 COMBINATIONS
                    
                    case ( F.B_N + F.B_E )  % NORTH EAST BOUNDARY
 
                        Cpo(j-1,i) = Cpo(j-1,i) * ( 1 + ((m_p/D_p)*(Uo(j-1,i)-sigma)) );                          
                        Cpo(j,i+1) = Cpo(j,i+1) * ( 1 + ((m_p/D_p)*(Uo(j,i+1)-sigma)) );
 
                        Cno(j-1,i) = Cno(j-1,i) * ( 1 - ((m_n/D_n)*(Uo(j-1,i)-sigma)) );                        
                        Cno(j,i+1) = Cno(j,i+1) * ( 1 - ((m_n/D_n)*(Uo(j,i+1)-sigma)) );                                                
                            
                    case ( F.B_N + F.B_W )  % NORTH WEST BOUNDARY
                        
                        Cpo(j-1,i) = Cpo(j-1,i) * ( 1 + ((m_p/D_p)*(Uo(j-1,i)-sigma)) );                         
                        Cpo(j,i-1) = Cpo(j,i-1) * ( 1 + ((m_p/D_p)*(Uo(j,i-1)-sigma)) );
  
                        Cno(j-1,i) = Cno(j-1,i) * ( 1 - ((m_n/D_n)*(Uo(j-1,i)-sigma)) );                         
                        Cno(j,i-1) = Cno(j,i-1) * ( 1 - ((m_n/D_n)*(Uo(j,i-1)-sigma)) ); 
                        
                    case ( F.B_S + F.B_E )  % SOUTH EAST BOUNDARY
                        
                        Cpo(j+1,i) = Cpo(j+1,i) * ( 1 + ((m_p/D_p)*(Uo(j+1,i)-sigma)) );                         
                        Cpo(j,i+1) = Cpo(j,i+1) * ( 1 + ((m_p/D_p)*(Uo(j,i+1)-sigma)) );
                         
                        Cno(j+1,i) = Cno(j+1,i) * ( 1 - ((m_n/D_n)*(Uo(j+1,i)-sigma)) );                          
                        Cno(j,i+1) = Cno(j,i+1) * ( 1 - ((m_n/D_n)*(Uo(j,i+1)-sigma)) );                                                
                          
                    case ( F.B_S + F.B_W )  % SOUTH WEST BOUNDARY
                        
                        Cpo(j+1,i) = Cpo(j+1,i) * ( 1 + ((m_p/D_p)*(Uo(j+1,i)-sigma)) );                           
                        Cpo(j,i-1) = Cpo(j,i-1) * ( 1 + ((m_p/D_p)*(Uo(j,i-1)-sigma)) );

                        Cno(j+1,i) = Cno(j+1,i) * ( 1 - ((m_n/D_n)*(Uo(j+1,i)-sigma)) );                         
                        Cno(j,i-1) = Cno(j,i-1) * ( 1 - ((m_n/D_n)*(Uo(j,i-1)-sigma)) );                                                
                                                    
                    case ( F.B_N + F.B_S )  % NORTH SOUTH BOUNDARY    
                      
                        Cpo(j-1,i) = Cpo(j-1,i) * ( 1 + ((m_p/D_p)*(Uo(j-1,i)-sigma)) );                        
                        Cpo(j+1,i) = Cpo(j+1,i) * ( 1 + ((m_p/D_p)*(Uo(j+1,i)-sigma)) );                           

                        Cno(j-1,i) = Cno(j-1,i) * ( 1 - ((m_n/D_n)*(Uo(j-1,i)-sigma)) ); 
                        Cno(j+1,i) = Cno(j+1,i) * ( 1 - ((m_n/D_n)*(Uo(j+1,i)-sigma)) );                         
                                               
                    case ( F.B_E + F.B_W )  % EAST WEST BOUNDARY    
                        
                        Cpo(j,i+1) = Cpo(j,i+1) * ( 1 + ((m_p/D_p)*(Uo(j,i+1)-sigma)) );                           
                        Cpo(j,i-1) = Cpo(j,i-1) * ( 1 + ((m_p/D_p)*(Uo(j,i-1)-sigma)) );

                        Cno(j,i+1) = Cno(j,i+1) * ( 1 - ((m_n/D_n)*(Uo(j,i+1)-sigma)) );                         
                        Cno(j,i-1) = Cno(j,i-1) * ( 1 - ((m_n/D_n)*(Uo(j,i-1)-sigma)) );                         
                        
                        
                    % 3 DIRECTIONS: TOTAL 4 COMBINATIONS   
                                                
                    case ( F.B_N + F.B_E + F.B_W )  % NORTH EAST WEST BOUNDARY
                        
                        Cpo(j-1,i) = Cpo(j-1,i) * ( 1 + ((m_p/D_p)*(Uo(j-1,i)-sigma)) );  
                        Cpo(j,i+1) = Cpo(j,i+1) * ( 1 + ((m_p/D_p)*(Uo(j,i+1)-sigma)) );
                        Cpo(j,i-1) = Cpo(j,i-1) * ( 1 + ((m_p/D_p)*(Uo(j,i-1)-sigma)) );                        
 
                        Cno(j-1,i) = Cno(j-1,i) * ( 1 - ((m_n/D_n)*(Uo(j-1,i)-sigma)) );                         
                        Cno(j,i+1) = Cno(j,i+1) * ( 1 - ((m_n/D_n)*(Uo(j,i+1)-sigma)) );                                                
                        Cno(j,i-1) = Cno(j,i-1) * ( 1 - ((m_p/D_p)*(Uo(j,i-1)-sigma)) );        
                        
                    case ( F.B_S + F.B_E + F.B_W )  % SOUTH EAST WEST BOUNDARY
                        
                        Cpo(j+1,i) = Cpo(j+1,i) * ( 1 + ((m_p/D_p)*(Uo(j+1,i)-sigma)) );  
                        Cpo(j,i+1) = Cpo(j,i+1) * ( 1 + ((m_p/D_p)*(Uo(j,i+1)-sigma)) );
                        Cpo(j,i-1) = Cpo(j,i-1) * ( 1 + ((m_p/D_p)*(Uo(j,i-1)-sigma)) );                        
 
                        Cno(j+1,i) = Cno(j+1,i) * ( 1 - ((m_n/D_n)*(Uo(j+1,i)-sigma)) );                         
                        Cno(j,i+1) = Cno(j,i+1) * ( 1 - ((m_n/D_n)*(Uo(j,i+1)-sigma)) );                                                
                        Cno(j,i-1) = Cno(j,i-1) * ( 1 - ((m_p/D_p)*(Uo(j,i-1)-sigma)) );                         
                                                    
                    case ( F.B_N + F.B_S + F.B_W )  % NORTH SOUTH WEST BOUNDARY
                        
                        Cpo(j-1,i) = Cpo(j-1,i) * ( 1 + ((m_p/D_p)*(Uo(j-1,i)-sigma)) );                        
                        Cpo(j+1,i) = Cpo(j+1,i) * ( 1 + ((m_p/D_p)*(Uo(j+1,i)-sigma)) );  
                        Cpo(j,i-1) = Cpo(j,i-1) * ( 1 + ((m_p/D_p)*(Uo(j,i-1)-sigma)) );                        
 
                        Cno(j-1,i) = Cno(j-1,i) * ( 1 - ((m_n/D_n)*(Uo(j-1,i)-sigma)) );                           
                        Cno(j+1,i) = Cno(j+1,i) * ( 1 - ((m_n/D_n)*(Uo(j+1,i)-sigma)) );                                                                     
                        Cno(j,i-1) = Cno(j,i-1) * ( 1 - ((m_p/D_p)*(Uo(j,i-1)-sigma)) );                          

                    case ( F.B_N + F.B_S + F.B_E )  % NORTH SOUTH EAST BOUNDARY

                        Cpo(j-1,i) = Cpo(j-1,i) * ( 1 + ((m_p/D_p)*(Uo(j-1,i)-sigma)) );                        
                        Cpo(j+1,i) = Cpo(j+1,i) * ( 1 + ((m_p/D_p)*(Uo(j+1,i)-sigma)) );  
                        Cpo(j,i+1) = Cpo(j,i+1) * ( 1 + ((m_p/D_p)*(Uo(j,i+1)-sigma)) );                        

                        Cno(j-1,i) = Cno(j-1,i) * ( 1 - ((m_n/D_n)*(Uo(j-1,i)-sigma)) );                           
                        Cno(j+1,i) = Cno(j+1,i) * ( 1 - ((m_n/D_n)*(Uo(j+1,i)-sigma)) );                                                                     
                        Cno(j,i+1) = Cno(j,i+1) * ( 1 - ((m_p/D_p)*(Uo(j,i+1)-sigma)) );                         
                        
                    otherwise
                        
                end % END OF SWITCH
                
            end % END OF IF
            
            
            if ( (F.LAYER(j,i) == F.OBJECT) & (P.BC_CONC == 2) )
    
                switch ( F.FLAG(j,i) )
                         
                    % 1 DIRECTIONS: TOTAL 4 COMBINATIONS
                    
                    case ( F.B_E )  % EAST BOUNDARY
                        
                        Cpo(j,i) = Cpo(j,i+1)*(1 + (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p))/(1 - (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p)); 
                        
                        Cno(j,i) = Cno(j,i+1)*(1 - (m_n*(Uo(j,i+1)-Uo(j,i)))/(2*D_n))/(1 + (m_n*(Uo(j,i+1)-Uo(j,i)))/(2*D_n));                         
     
                    case ( F.B_W )  % WEST BOUNDARY
                        
                        Cpo(j,i) = Cpo(j,i-1)*(1 - (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p)); 
                        
                        Cno(j,i) = Cno(j,i-1)*(1 + (m_n*(Uo(j,i)-Uo(j,i-1)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j,i-1)))/(2*D_n));                         
                                               
                    case ( F.B_N )  % NORTH BOUNDARY
                        
                        Cpo(j,i) = Cpo(j-1,i)*(1 - (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p)); 
                        
                        Cno(j,i) = Cno(j-1,i)*(1 + (m_n*(Uo(j,i)-Uo(j-1,i)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j-1,i)))/(2*D_n));                         
                
                    case ( F.B_S )  % SOUTH BOUNDARY
                         
                        Cpo(j,i) = Cpo(j+1,i)*(1 - (m_p*(Uo(j,i)-Uo(j+1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j+1,i)))/(2*D_p)); 
                        
                        Cno(j,i) = Cno(j+1,i)*(1 + (m_n*(Uo(j,i)-Uo(j+1,i)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j+1,i)))/(2*D_n)); 
                        
                        
                    % 2 DIRECTIONS: TOTAL 6 COMBINATIONS
                    
                    case ( F.B_N + F.B_E )  % NORTH EAST BOUNDARY
 
                        Cpo(j,i) = (  Cpo(j-1,i)*(1 - (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p)) ...
                                    + Cpo(j,i+1)*(1 - (m_p*(Uo(j,i)-Uo(j,i+1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i+1)))/(2*D_p)) )/2; 
                        
                        Cno(j,i) = (  Cno(j-1,i)*(1 + (m_n*(Uo(j,i)-Uo(j-1,i)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j-1,i)))/(2*D_n)) ...
                                    + Cno(j,i+1)*(1 + (m_n*(Uo(j,i)-Uo(j,i+1)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j,i+1)))/(2*D_n)) )/2; 
                            
                    case ( F.B_N + F.B_W )  % NORTH WEST BOUNDARY
 
                        Cpo(j,i) = (  Cpo(j-1,i)*(1 - (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p)) ...
                                    + Cpo(j,i-1)*(1 - (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p)) )/2; 
                        
                        Cno(j,i) = (  Cno(j-1,i)*(1 + (m_n*(Uo(j,i)-Uo(j-1,i)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j-1,i)))/(2*D_n)) ...
                                    + Cno(j,i-1)*(1 + (m_n*(Uo(j,i)-Uo(j,i-1)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j,i-1)))/(2*D_n)) )/2;                        
                     
                    case ( F.B_S + F.B_E )  % SOUTH EAST BOUNDARY

                        Cpo(j,i) = (  Cpo(j+1,i)*(1 - (m_p*(Uo(j,i)-Uo(j+1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j+1,i)))/(2*D_p)) ...
                                    + Cpo(j,i+1)*(1 - (m_p*(Uo(j,i)-Uo(j,i+1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i+1)))/(2*D_p)) )/2; 
                        
                        Cno(j,i) = (  Cno(j+1,i)*(1 + (m_n*(Uo(j,i)-Uo(j+1,i)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j+1,i)))/(2*D_n)) ...
                                    + Cno(j,i+1)*(1 + (m_n*(Uo(j,i)-Uo(j,i+1)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j,i+1)))/(2*D_n)) )/2;                         
                                                                         
                    case ( F.B_S + F.B_W )  % SOUTH WEST BOUNDARY

                        Cpo(j,i) = (  Cpo(j+1,i)*(1 - (m_p*(Uo(j,i)-Uo(j+1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j+1,i)))/(2*D_p)) ...
                                    + Cpo(j,i-1)*(1 - (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p)) )/2; 
                        
                        Cno(j,i) = (  Cno(j+1,i)*(1 + (m_n*(Uo(j,i)-Uo(j+1,i)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j+1,i)))/(2*D_n)) ...
                                    + Cno(j,i-1)*(1 + (m_n*(Uo(j,i)-Uo(j,i-1)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j,i-1)))/(2*D_n)) )/2; 
                                                                                                 
                    case ( F.B_N + F.B_S )  % NORTH SOUTH BOUNDARY    
                        
                        Cpo(j,i) = (  Cpo(j-1,i)*(1 - (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p)) ... 
                                    + Cpo(j+1,i)*(1 - (m_p*(Uo(j,i)-Uo(j+1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j+1,i)))/(2*D_p)) )/2; 
                        
                        Cno(j,i) = (  Cno(j-1,i)*(1 + (m_n*(Uo(j,i)-Uo(j-1,i)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j-1,i)))/(2*D_n)) ...                       
                                    + Cno(j+1,i)*(1 + (m_n*(Uo(j,i)-Uo(j+1,i)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j+1,i)))/(2*D_n)) )/2;                       
                                                                      
                    case ( F.B_E + F.B_W )  % EAST WEST BOUNDARY    
                        
                        Cpo(j,i) = (  Cpo(j,i+1)*(1 - (m_p*(Uo(j,i)-Uo(j,i+1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i+1)))/(2*D_p)) ... 
                                    + Cpo(j,i-1)*(1 - (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p)) )/2; 
                        
                        Cno(j,i) = (  Cno(j,i-1)*(1 + (m_n*(Uo(j,i)-Uo(j,i-1)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j,i-1)))/(2*D_n)) ...                         
                                   +  Cno(j,i+1)*(1 + (m_n*(Uo(j,i)-Uo(j,i+1)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j,i+1)))/(2*D_n)) )/2;                         
                        
                        
                    % 3 DIRECTIONS: TOTAL 4 COMBINATIONS   
                                                
                    case ( F.B_N + F.B_E + F.B_W )  % NORTH EAST WEST BOUNDARY

                        Cpo(j,i) = (  Cpo(j-1,i)*(1 - (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p)) ...
                                    + Cpo(j,i+1)*(1 - (m_p*(Uo(j,i)-Uo(j,i+1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i+1)))/(2*D_p)) ...
                                    + Cpo(j,i-1)*(1 - (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p)) )/3; 
                        
                        Cno(j,i) = (  Cno(j-1,i)*(1 + (m_n*(Uo(j,i)-Uo(j-1,i)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j-1,i)))/(2*D_n)) ...
                                    + Cno(j,i+1)*(1 + (m_n*(Uo(j,i)-Uo(j,i+1)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j,i+1)))/(2*D_n)) ...
                                    + Cno(j,i-1)*(1 + (m_n*(Uo(j,i)-Uo(j,i-1)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j,i-1)))/(2*D_n)) )/3; 
     
                    case ( F.B_S + F.B_E + F.B_W )  % SOUTH EAST WEST BOUNDARY
                        
                        Cpo(j,i) = (  Cpo(j+1,i)*(1 - (m_p*(Uo(j,i)-Uo(j+1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j+1,i)))/(2*D_p)) ...
                                    + Cpo(j,i+1)*(1 - (m_p*(Uo(j,i)-Uo(j,i+1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i+1)))/(2*D_p)) ...
                                    + Cpo(j,i-1)*(1 - (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p)) )/3; 
                        
                        Cno(j,i) = (  Cno(j+1,i)*(1 + (m_n*(Uo(j,i)-Uo(j+1,i)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j+1,i)))/(2*D_n)) ...
                                    + Cno(j,i+1)*(1 + (m_n*(Uo(j,i)-Uo(j,i+1)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j,i+1)))/(2*D_n)) ...
                                    + Cno(j,i-1)*(1 + (m_n*(Uo(j,i)-Uo(j,i-1)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j,i-1)))/(2*D_n)) )/3;                        
                                                                           
                    case ( F.B_N + F.B_S + F.B_W )  % NORTH SOUTH WEST BOUNDARY

                        Cpo(j,i) = (  Cpo(j-1,i)*(1 - (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p)) ...
                                    + Cpo(j+1,i)*(1 - (m_p*(Uo(j,i)-Uo(j+1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j+1,i)))/(2*D_p)) ...
                                    + Cpo(j,i-1)*(1 - (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p)) )/3; 
                                
                        Cno(j,i) = (  Cno(j-1,i)*(1 + (m_n*(Uo(j,i)-Uo(j-1,i)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j-1,i)))/(2*D_n)) ...
                                    + Cno(j+1,i)*(1 + (m_n*(Uo(j,i)-Uo(j+1,i)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j+1,i)))/(2*D_n)) ...
                                    + Cno(j,i-1)*(1 + (m_n*(Uo(j,i)-Uo(j,i-1)))/(2*D_n))/(1 - (m_n*(Uo(j,i)-Uo(j,i-1)))/(2*D_n)) )/3;                                 
                        
                    case ( F.B_N + F.B_S + F.B_E )  % NORTH SOUTH EAST BOUNDARY

                        Cpo(j,i) = (  Cpo(j-1,i)*(1 - (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p)) ...
                                    + Cpo(j+1,i)*(1 - (m_p*(Uo(j,i)-Uo(j+1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j+1,i)))/(2*D_p)) ...
                                    + Cpo(j,i+1)*(1 - (m_p*(Uo(j,i)-Uo(j,i+1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i+1)))/(2*D_p)) )/3; 

                        Cno(j,i) = (  Cno(j-1,i)*(1 + (m_n*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))/(1 - (m_n*(Uo(j,i)-Uo(j-1,i)))/(2*D_n)) ...
                                    + Cno(j+1,i)*(1 + (m_n*(Uo(j,i)-Uo(j+1,i)))/(2*D_p))/(1 - (m_n*(Uo(j,i)-Uo(j+1,i)))/(2*D_n)) ...
                                    + Cno(j,i+1)*(1 + (m_n*(Uo(j,i)-Uo(j,i+1)))/(2*D_p))/(1 - (m_n*(Uo(j,i)-Uo(j,i+1)))/(2*D_n)) )/3;                                 
                                
                    otherwise
                        
                end % END OF SWITCH
                
            end % END OF IF            
                        
        end % END OF FOR

    end % END OF FOR
  
end
