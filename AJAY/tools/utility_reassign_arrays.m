
message  = 'Enter the last iteration? ';
timestep = input(message)

% REASSIGN TIMESTEPS
iteration_t = timestep;
i = timestep;

% residual_U_temp  = P.residual_U(1:iteration_t); 
% iteration_U_temp = P.iteration_U(1:iteration_t);
% 
% P.residual_U  = residual_U_temp;
% P.iteration_U = iteration_U_temp;

clear residual_U_temp iteration_U_temp

% U, C_p and C_n ARRAYS
U_temp  = Ut(:,:,[1:timestep]);
Cp_temp = Cpt(:,:,[1:timestep]);
Cn_temp = Cnt(:,:,[1:timestep]);

Ut  = U_temp;
Cpt = Cp_temp;
Cnt = Cn_temp;

clear U_temp Cp_temp Cn_temp