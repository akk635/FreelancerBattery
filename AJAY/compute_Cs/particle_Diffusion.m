function [Csn] = particle_Diffusion(Csn, Cso, P, F )   

    % ficks law of dissusion in both X, Y directions
    for i = 2 : P.Nx+1
        for j = 2: P.Ny+1
            if (F.FLAG(j,i) == 0)
                Csn(j,i) = Cso(j,i) + P.D_s * P.dt * (  ((Cso(j,i+1) - 2*Cso(j,i) + Cso(j,i-1))/(P.hx * P.hx)) ...
                                        + ((Cso(j+1,i) - 2*Cso(j,i) + Cso(j-1,i))/(P.hy * P.hy)) );
            end
            if (F.FLAG(j,i) > 0 & F.FLAG(j,i) < 16)
                Csn(j,i) = Cso(j,i);
            end
        end
    end
   
end