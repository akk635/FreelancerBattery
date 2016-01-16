function [Csn] = particle_Diffusion(P, F, t, Csn, Cso, Cpn, i_app )
    
    % scale C_p
    realCpn = Cpn /P.NA/P.Debye^3;
    
    % i_rn = k * F * C_p^0.5

end