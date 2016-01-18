function [Cso] = particle_Boundary(Cso, i_app, P, F, Cpn, C_smax)
    
    % to determine the reaction constant
    % hypothesis i_rxn = k_rxn * F * Cpn^0.5 * Cso^0.5 * (C_smax - Cso)^0.5
      
    totCurr = 0;  % total chrg transfer current
    % boundary representation F_N = 1, F_S = 2, F_E = 4, F_W = 8   
    %initially lets only consider the intercalation of positive ions
    for i = 2 : P.Nx+1
        for j = 2: P.Ny+1
            if ((F.LAYER(j, i) == F.OBJECT) & (F.FLAG(j,i) > 0))
                % determine the direction of neighboring fluid cell
                bndBits = dec2bin(F.FLAG(j,i),4);

                if (bndBits(4) == 1) % F_N
                    totCurr = totCurr + (F * Cso(j,i)^0.5 * (Cpn(j-1,i)/P.NA/P.Debye)^0.5 * (C_smax - Cso(j,i))^0.5)*P.hx;
                end
                if (bndBits(3) == 1) % F_S
                    totCurr = totCurr + (F * Cso(j,i)^0.5 * (Cpn(j+1,i)/P.NA/P.Debye)^0.5 * (C_smax - Cso(j,i))^0.5)*P.hx;
                end
                if (bndBits(2) == 1) % F_E
                    totCurr = totCurr + (F * Cso(j,i)^0.5 * (Cpn(j,i+1)/P.NA/P.Debye)^0.5 * (C_smax - Cso(j,i))^0.5)*P.hy;
                end
                if (bndBits(1) == 1) % F_W
                    totCurr = totCurr + (F * Cso(j,i)^0.5 * (Cpn(j,i-1)/P.NA/P.Debye)^0.5 * (C_smax - Cso(j,i))^0.5)*P.hy;
                end               
            end
        end
    end
    
    % now the calculation of k_rxn
    k_rxn = i_app * P.Ly/totCurr;
    % then the actual boundary current
    for i = 2 : P.Nx+1
        for j = 2: P.Ny+1
            
            if ((F.LAYER(j, i) == F.OBJECT) & (F.FLAG(j,i) > 0))
                % determine the direction of neighboring fluid cell
                bndBits = dec2bin(F.FLAG(j,i),4);

                if (bndBits(4) == 1) % F_N
                    Cso(j,i) = Cso(j+1,i) + (k_rxn * Cso(j,i)^0.5 * (Cpn(j-1,i)/P.NA/P.Debye)^0.5 * (C_smax - Cso(j,i))^0.5 * P.hy / P.D_s);
                end
                if (bndBits(3) == 1) % F_S
                    Cso(j,i) = Cso(j-1,i) + (k_rxn * Cso(j,i)^0.5 * (Cpn(j+1,i)/P.NA/P.Debye)^0.5 * (C_smax - Cso(j,i))^0.5 * P.hy / P.D_s);
                end
                if (bndBits(2) == 1) % F_E
                    Cso(j,i) = Cso(j,i-1) + (k_rxn * Cso(j,i)^0.5 * (Cpn(j,i+1)/P.NA/P.Debye)^0.5 * (C_smax - Cso(j,i))^0.5 * P.hx / P.D_s);                    
                end
                if (bndBits(1) == 1) % F_W
                    Cso(j,i) = Cso(j,i+1) + (k_rxn * Cso(j,i)^0.5 * (Cpn(j,i-1)/P.NA/P.Debye)^0.5 * (C_smax - Cso(j,i))^0.5 * P.hx / P.D_s);
                end               
            end
        end
    end
    
end