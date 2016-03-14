
function [ P S D F Uo Un Ut Cpo Cpn Cpt Cno Cnn Cnt Cso Csn Cst] = init( filename )
  
  % P: A STRUCTURE FOR MODEL PARAMETERS
  [ P S D ] = parameters( filename );

  % F: A STRUCTURE FOR FLAG FIELD
  [ F ] = init_flag( filename, P );  
  
  % TODO: BETTER WAY TO STORE TENSOR
  [ Uo Un Ut Cpo Cpn Cpt Cno Cnn Cnt Cso Csn Cst] = init_matrices( P, S, D, F );

end
