
function [ P ] = utility_store_residual_iteration_U( P, timestep, iteration, residual )

    P.residual_U(timestep)  = residual;
    P.iteration_U(timestep) = iteration;

end 