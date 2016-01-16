
function [ F ] = init_flag( filename, P )

    % FLUID & BOUNDARY CELLS
    % USE: TO DETERMINE WHICH KIND OF CELLS SURROND A GIVEN CELL 
    F.C_F  = 16; % FLUID CELL
    F.C_B  = 0;  % BOUNDARY CELL

    % Fluid boundaries
    F.B_N  = 1; % BOUNDARY NORTH
    F.B_S  = 2; % BOUNDARY SOUTH 
    F.B_E  = 4; % BOUNDARY EAST
    F.B_W  = 8; % BOUNDARY WEST
    
    % OBJECT CELL, DIFFUSE LAYER CELL & STERN LAYER CELL
    F.OBJECT  = 0;  % OBJECT CELL
    F.DIFFUSE = 1;  % DIFFUSE LAYER CELL
    F.STERN   = 2;  % STERN LAYER CELL 
    F.THIRD   = 3;  % THICK STERN LAYER CELL
    F.FOURTH  = 4;  % THICKER STERN LAYER CELL 
    
    % LOADS THE IMAGE I.E. BITMAP (0S & 1S)
    S = load( filename );
    F.IMAGE = sparse(S.A);  clear S;
    
    [F.Nx,F.Ny] = size(F.IMAGE);    
    
    F.GEO         = sparse(F.Nx+2,F.Ny+2);    % GEOMETRY MATRIX
    F.FLAG        = sparse(F.Nx+2,F.Ny+2);    % FLAG MATRIX    
    F.LAYER       = sparse(F.Nx+2,F.Ny+2);    % LAYER MATRIX  
    F.LAYER_CLONE = sparse(F.Nx+2,F.Ny+2);    
    
    % COMPUTE GEOMETRY MATRIX
    F.GEO([2:F.Nx+1],[2:F.Ny+1]) = F.IMAGE;  
    
    % COMPUTE FLAG MATRIX
    F.FLAG           = F.C_F*F.GEO;   
    
    F.FLAG(1,:)      = F.B_N;
    F.FLAG(F.Nx+2,:) = F.B_S;
    F.FLAG(:,1)      = F.B_W;
    F.FLAG(:,F.Ny+2) = F.B_E;
    
    % allocating boundary cdns to the fluid cells
    for j = 2:F.Ny+1
        
        for i = 2:F.Nx+1
            
            if ( F.FLAG(j+1,i) >= F.C_F )
                F.FLAG(j,i) = F.FLAG(j,i) + F.B_S;                   
            end
            
            if ( F.FLAG(j-1,i) >= F.C_F )
                F.FLAG(j,i) = F.FLAG(j,i) + F.B_N;                
            end

            if ( F.FLAG(j,i+1) >= F.C_F )
                F.FLAG(j,i) = F.FLAG(j,i) + F.B_E;                    
            end

            if ( F.FLAG(j,i-1) >= F.C_F )
                F.FLAG(j,i) = F.FLAG(j,i) + F.B_W;                   
            end
                
        end
            
    end
    
    % COMPUTE LAYER MATRIX
    F.LAYER = F.GEO;
    
    if ( P.STERN == true )
    
%         % STERN LAYER
%         for j = 2:F.Ny+1
%         
%             for i = 2:F.Nx+1
%                 
%                 if ( F.FLAG(j,i) < F.C_F )
%                     
%                     if ( F.FLAG(j,i+1) > F.C_F )
%                         F.LAYER(j,i+1) = F.STERN;
%                     end
% 
%                     if ( F.FLAG(j,i-1) > F.C_F )
%                         F.LAYER(j,i-1) = F.STERN;
%                     end                    
%                     
%                     if ( F.FLAG(j+1,i) > F.C_F )
%                         F.LAYER(j+1,i) = F.STERN;
%                     end                     
%                     
%                     if ( F.FLAG(j-1,i) > F.C_F )
%                         F.LAYER(j-1,i) = F.STERN;
%                     end   
%                     
%                     if ( F.FLAG(j+1,i+1) > F.C_F )
%                         F.LAYER(j+1,i+1) = F.STERN;                    
%                     end
% 
%                     if ( F.FLAG(j-1,i+1) > F.C_F )
%                         F.LAYER(j-1,i+1) = F.STERN;                    
%                     end
%                     
%                     if ( F.FLAG(j+1,i-1) > F.C_F )
%                         F.LAYER(j+1,i-1) = F.STERN;                   
%                     end                    
%                     
%                     if ( F.FLAG(j-1,i-1) > F.C_F )
%                         F.LAYER(j-1,i-1) = F.STERN;                    
%                     end
%                 
%                 end
%                 
%             end
%             
%         end
        
%         % THICK STERN LAYER   
%         for j = 2:F.Ny+1
%         
%             for i = 2:F.Nx+1
%                 
%                 if ( F.LAYER(j,i) == F.STERN )
%                     
%                     if ( (F.FLAG(j,i+1) > F.C_F) & (F.LAYER(j,i+1) ~= F.STERN) ) 
%                         F.LAYER(j,i+1) = F.THIRD;
%                     end
% 
%                     if ( (F.FLAG(j,i-1) > F.C_F) & (F.LAYER(j,i-1) ~= F.STERN) )
%                         F.LAYER(j,i-1) = F.THIRD;
%                     end                    
%                     
%                     if ( (F.FLAG(j+1,i) > F.C_F) & (F.LAYER(j+1,i) ~= F.STERN) )
%                         F.LAYER(j+1,i) = F.THIRD;
%                     end                     
%                     
%                     if ( (F.FLAG(j-1,i) > F.C_F) & (F.LAYER(j-1,i) ~= F.STERN) )
%                         F.LAYER(j-1,i) = F.THIRD;
%                     end   
%                     
%                     if ( (F.FLAG(j+1,i+1) > F.C_F) & (F.LAYER(j+1,i+1) ~= F.STERN) )
%                         F.LAYER(j+1,i+1) = F.THIRD;                    
%                     end
% 
%                     if ( (F.FLAG(j-1,i+1) > F.C_F) & (F.LAYER(j-1,i+1) ~= F.STERN) )
%                         F.LAYER(j-1,i+1) = F.THIRD;                    
%                     end
%                     
%                     if ( (F.FLAG(j+1,i-1) > F.C_F) & (F.LAYER(j+1,i-1) ~= F.STERN) )
%                         F.LAYER(j+1,i-1) = F.THIRD;                   
%                     end                    
%                     
%                     if ( (F.FLAG(j-1,i-1) > F.C_F) & (F.LAYER(j-1,i-1) ~= F.STERN) )
%                         F.LAYER(j-1,i-1) = F.THIRD;                    
%                     end
%                 
%                 end
%                 
%             end
%             
%         end 
        
        F.LAYER_CLONE = F.LAYER;
        
    end
    
end
