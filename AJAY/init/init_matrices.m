
function [ Uo Un Ut Cpo Cpn Cpt Cno Cnn Cnt ] = init_matrices( P, S, D, F )

  % GENERATING INITIAL CONCENTRATION CP, CN AND U
  Uo  = sparse(P.Nx+2,P.Ny+2);
  Un  = sparse(P.Nx+2,P.Ny+2);
  Un(:,:) = 0;
  Uo(:,:) = 0;  
  % NOTE: 
  % 1. DIFFUSE LAYER IS SOLVED FIRST
  % 2. THE RESULTS FROM DIFFUSE LAYER ARE TAKEN AS THE INITIAL CONDITIONS
  %    FOR CP & CN IN THE STERN LAYER
  % 3. THUS, NO NEED TO ASSIGN CP & CN FOR THE STERN LAYER
  Cpo = D.Co*F.GEO; 
  Cpn = D.Co*F.GEO; 
  
  Cno = D.Co*F.GEO;
  Cnn = D.Co*F.GEO;

  Ut  = zeros(P.Nx+2,P.Ny+2,200);
  Cpt = zeros(P.Nx+2,P.Ny+2,200);
  Cnt = zeros(P.Nx+2,P.Ny+2,200);
  
  Csn = sparse(P.Nx+2, P.Ny+2);
  Cso = sparse(P.Nx+2, P.Ny+2);
  Csn(:,:) = 0;
  Cso(:,:) = 0;
  Cso(F.GEO == 0) = 0.2 * P.C_smax;
  
  Cst = zeros(P.Nx+2,P.Ny+2,200);
end
