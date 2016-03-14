function [Cso, Cpo] = particle_Boundary(Cso, Cpo, Uo, j_rxn, P, F, D)
    
    % to determine the reaction constant
    % hypothesis i_rxn = k_rxn * F * Cpn^0.5 * Cso^0.5 * (C_smax - Cso)^0.5
    % in this case considering the current to be going inside the particle
      
%   totCurr = 0;  % total chrg transfer current
    C_smax = P.C_smax;
    k_rxn = P.k_rxn;
%   ND equations scaling factors    
    D_p   = D.D_p;
    D_n   = D.D_n;
    m_p   = D.m_p;
    m_n   = D.m_n;
  
    Nx = P.Nx;
    Ny = P.Ny;
    hx = D.hx;
    hy = D.hy;
    dt = D.dt;
    
%     boundary representation F_N = 1, F_S = 2, F_E = 4, F_W = 8   
%     initially lets only consider the intercalation of positive ions
%     for i = 2 : P.Nx+1
%         for j = 2: P.Ny+1
%             if ((F.LAYER(j, i) == F.OBJECT) & (F.FLAG(j,i) > 0))
%                 determine the direction of neighboring fluid cell
%                 bndBits = dec2bin(F.FLAG(j,i),4);
% 
%                 if (bndBits(4) == 1) % F_N
%                     totCurr = totCurr + (F * Cso(j,i)^0.5 * (Cpn(j-1,i)/P.NA/P.Debye^3)^0.5 * (C_smax - Cso(j,i))^0.5)*P.hx;
%                 end
%                 if (bndBits(3) == 1) % F_S
%                     totCurr = totCurr + (F * Cso(j,i)^0.5 * (Cpn(j+1,i)/P.NA/P.Debye^3)^0.5 * (C_smax - Cso(j,i))^0.5)*P.hx;
%                 end
%                 if (bndBits(2) == 1) % F_E
%                     totCurr = totCurr + (F * Cso(j,i)^0.5 * (Cpn(j,i+1)/P.NA/P.Debye^3)^0.5 * (C_smax - Cso(j,i))^0.5)*P.hy;
%                 end
%                 if (bndBits(1) == 1) % F_W
%                     totCurr = totCurr + (F * Cso(j,i)^0.5 * (Cpn(j,i-1)/P.NA/P.Debye^3)^0.5 * (C_smax - Cso(j,i))^0.5)*P.hy;
%                 end               
%             end
%         end
%     end
%     
%     now the calculation of k_rxn
%     k_rxn = i_app * P.Ly/totCurr;

%   flux at the electrode-electrolyte interface j_rxn/F
%   ND Equation Dp * \nabla C_p + mp C_p \nabla  U = j_rxn * Dp/F
    % then the actual boundary current
    tempCpW = 0.0;
    tempCpE = 0.0;
    for i = 2 : P.Nx+1
        for j = 2: P.Ny+1
            
            if ((F.LAYER(j, i) == F.OBJECT) & (F.FLAG(j,i) > 0))
               
                % determine the direction of neighboring fluid cell
                bndBits = dec2bin(F.FLAG(j,i),4);                

                if (bndBits(4) == '1') % F_N
                    j_rxn(1) = k_rxn * Cso(j,i)^0.5 * (Cpo(j-1,i)/P.NA/P.Debye^3)^0.5 * (C_smax - Cso(j,i))^0.5;                   
                end
                if (bndBits(3) == '1') % F_S
                    j_rxn(2) = k_rxn * Cso(j,i)^0.5 * ( Cpo(j+1,i)/P.NA/P.Debye^3 )^0.5 * (C_smax - Cso(j,i))^0.5;                    
                end
                if (bndBits(2) == '1') % F_E
                    j_rxn(3) = k_rxn * Cso(j,i)^0.5 * (Cpo(j,i+1)/P.NA/P.Debye^3)^0.5 * (C_smax - Cso(j,i))^0.5;
                    if (j_rxn(3) > tempCpE)
                        tempCpE = j_rxn(3);
                    end
                end
                if (bndBits(1) == '1') % F_W
                    j_rxn(4) = k_rxn * Cso(j,i)^0.5 * (Cpo(j,i-1)/P.NA/P.Debye^3)^0.5 * (C_smax - Cso(j,i))^0.5;
                    if (j_rxn(4) > tempCpW)
                        tempCpW = j_rxn(4);
                    end
                end
                
                switch( F.FLAG(j,i) )
                    case (F.B_N)
                        Cso(j,i) = Cso(j+1,i) + (j_rxn(1) * P.hy /P.Faraday / P.D_s);
                        Cpo(j,i) = Cpo(j-1,i)*(1 - (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p)) ...
                                 - j_rxn(1) * hy/P.Faraday/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))* P.Debye^4 * P.NA /P.D;
                    case (F.B_S)
                        Cso(j,i) = Cso(j-1,i) + (j_rxn(2) * P.hy /P.Faraday / P.D_s);
                        Cpo(j,i) = Cpo(j+1,i)*(1 + (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p))/(1 - (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p)) ...
                                 - j_rxn(2) * hy/P.Faraday/(1 - (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p))* P.Debye^4 * P.NA /P.D;                        
                    case (F.B_E)
                        Cso(j,i) = Cso(j,i-1) + ( j_rxn(3) * P.hx /P.Faraday / P.D_s);
                        Cpo(j,i) = Cpo(j,i+1)*(1 + (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p))/(1 - (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p)) ...
                                   - j_rxn(3) * hx/P.Faraday/(1 - (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p))* P.Debye^4 * P.NA /P.D;
                    case (F.B_W)
                        Cso(j,i) = Cso(j,i+1) + (j_rxn(4) * P.hx / P.Faraday / P.D_s);
                        Cpo(j,i) = Cpo(j,i-1)*(1 - (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p)) ...
                                   - j_rxn(4) * hx / P.Faraday / (1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p))* P.Debye^4 * P.NA /P.D;   
                               
                    case (F.B_N + F.B_E)
                        Cso(j,i) = (Cso(j+1,i) + (j_rxn(1) * P.hy /P.Faraday / P.D_s) + Cso(j,i-1) + ( j_rxn(3) * P.hx /P.Faraday / P.D_s)) / 2;
                        Cpo(j,i) = ((Cpo(j-1,i)*(1 - (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p)) ...
                                    - j_rxn(1) * hy/P.Faraday/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))* P.Debye^4 * P.NA /P.D) + ...
                                    (Cpo(j,i+1)*(1 + (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p))/(1 - (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p)) ...
                                   - j_rxn(3) * hx/P.Faraday/(1 - (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p))* P.Debye^4 * P.NA /P.D))/2;
                    case (F.B_N + F.B_W)
                        Cso(j,i) = (Cso(j+1,i) + (j_rxn(1) * P.hy /P.Faraday / P.D_s) + Cso(j,i+1) + (j_rxn(4) * P.hx / P.Faraday / P.D_s))/2;
                        Cpo(j,i) = ((Cpo(j-1,i)*(1 - (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p)) ...
                                 - j_rxn(1) * hy/P.Faraday/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))* P.Debye^4 * P.NA /P.D) + ...
                                 (Cpo(j,i-1)*(1 - (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p)) ...
                                   - j_rxn(4) * hx / P.Faraday / (1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p))* P.Debye^4 * P.NA /P.D))/2;
                    case (F.B_S + F.B_E)
                        Cso(j,i) = (Cso(j-1,i) + (j_rxn(2) * P.hy /P.Faraday / P.D_s) + Cso(j,i-1) + ( j_rxn(3) * P.hx /P.Faraday / P.D_s))/2;
                        Cpo(j,i) = ((Cpo(j+1,i)*(1 + (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p))/(1 - (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p)) ...
                                 - j_rxn(2) * hy/P.Faraday/(1 - (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p))* P.Debye^4 * P.NA /P.D) + ...
                                 (Cpo(j,i+1)*(1 + (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p))/(1 - (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p)) ...
                                   - j_rxn(3) * hx/P.Faraday/(1 - (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p))* P.Debye^4 * P.NA /P.D))/2;
                    case (F.B_S + F.B_W)
                        Cso(j,i) = (Cso(j-1,i) + (j_rxn(2) * P.hy /P.Faraday / P.D_s) + Cso(j,i+1) + (j_rxn(4) * P.hx / P.Faraday / P.D_s))/2;
                        Cpo(j,i) = ((Cpo(j+1,i)*(1 + (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p))/(1 - (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p)) ...
                                 - j_rxn(2) * hy/P.Faraday/(1 - (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p))* P.Debye^4 * P.NA /P.D) + ...
                                 (Cpo(j,i-1)*(1 - (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p)) ...
                                   - j_rxn(4) * hx / P.Faraday / (1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p))* P.Debye^4 * P.NA /P.D))/2;
%                     case (F.B_N + F.B_E + F.B_S)
%                         Cso(j,i) = (Cso(j+1,i) + (j_rxn(1) * P.hy / P.D_s) + Cso(j,i-1) + ( j_rxn(3)* P.hx / P.D_s) + Cso(j-1,i) + (j_rxn(2) * P.hy / P.D_s))/3;
%                         Cpo(j,i) = (Cpo(j-1,i)*(1 - (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p)) ...
%                                     - j_rxn(1) * hy/P.Faraday/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p)) + ...
%                                     Cpo(j,i+1)*(1 + (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p))/(1 - (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p)) ...
%                                    - j_rxn(3) * hx/P.Faraday/(1 - (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p)) + ...
%                                    Cpo(j+1,i)*(1 + (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p))/(1 - (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p)) ...
%                                  - j_rxn(2) * hy/P.Faraday/(1 - (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p)))/3;
%                         
%                     case (F.B_N + F.B_W + F.B_S)
%                         Cso(j,i) = (Cso(j+1,i) + (j_rxn(1) * P.hy / P.D_s) + Cso(j,i+1) + (j_rxn(4) * P.hx / P.D_s) + Cso(j-1,i) + (j_rxn(2) * P.hy / P.D_s))/3;
%                         Cpo(j,i) = (Cpo(j-1,i)*(1 - (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p)) ...
%                                  - j_rxn(1) * hy/P.Faraday/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p)) + ...
%                                  Cpo(j,i-1)*(1 - (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p)) ...
%                                    - j_rxn(4) * hx / P.Faraday / (1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p)) + ...
%                                  Cpo(j+1,i)*(1 + (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p))/(1 - (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p)) ...
%                                    - j_rxn(2) * hy/P.Faraday/(1 - (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p)))/3;
%                     case (F.B_E + F.B_N + F.B_W)
%                         Cso(j,i) = (Cso(j+1,i) + (j_rxn(1) * P.hy / P.D_s) + Cso(j,i-1) + ( j_rxn(3)* P.hx / P.D_s) + Cso(j,i+1) + (j_rxn(4) * P.hx / P.D_s))/3 ;
%                         Cpo(j,i) = (Cpo(j-1,i)*(1 - (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p)) ...
%                                     - j_rxn(1) * hy/P.Faraday/(1 + (m_p*(Uo(j,i)-Uo(j-1,i)))/(2*D_p)) + ...
%                                     Cpo(j,i+1)*(1 + (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p))/(1 - (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p)) ...
%                                    - j_rxn(3) * hx/P.Faraday/(1 - (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p)) + ...
%                                     Cpo(j,i-1)*(1 - (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p)) ...
%                                    - j_rxn(4) * hx / P.Faraday / (1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p)))/3;
%                     case (F.B_E + F.B_S + F.B_W)
%                         Cso(j,i) = (Cso(j-1,i) + (j_rxn(2) * P.hy / P.D_s) + Cso(j,i-1) + ( j_rxn(3)* P.hx / P.D_s) + Cso(j,i+1) + (j_rxn(4) * P.hx / P.D_s))/3;                        
%                         Cpo(j,i) = (Cpo(j+1,i)*(1 + (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p))/(1 - (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p)) ...
%                                  - j_rxn(2) * hy/P.Faraday/(1 - (m_p*(Uo(j+1,i)-Uo(j,i)))/(2*D_p)) + ...
%                                  Cpo(j,i+1)*(1 + (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p))/(1 - (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p)) ...
%                                    - j_rxn(3) * hx/P.Faraday/(1 - (m_p*(Uo(j,i+1)-Uo(j,i)))/(2*D_p)) + ...
%                                    Cpo(j,i-1)*(1 - (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p))/(1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p)) ...
%                                    - j_rxn(4) * hx / P.Faraday / (1 + (m_p*(Uo(j,i)-Uo(j,i-1)))/(2*D_p)))/3;
                    otherwise
                end
                
            end
        end
    end    
end