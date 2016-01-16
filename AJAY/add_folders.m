
root   = pwd;
folder = {
          '',...
          'boundary_conditions'         , ...
          'boundary_conditions/diffuse' , ...
          'compute_CP_CN'               , ...
          'compute_CP_CN/diffuse'       , ...
          'compute_U'                   , ...
          'compute_U/diffuse'           , ...
          'EDL'                         , ...
          'EDL/diffuse'                 , ...
          'init'                        , ...
          'tools'                       , ...          
          '_images'                
         };

% ADDING FOLDER TO PATH
for i = 1:length(folder)
    path = strcat(root,'/',folder{i});
    addpath(path);
end

clear root folder path i
