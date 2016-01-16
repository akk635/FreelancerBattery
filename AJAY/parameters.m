
function [ P S D ] = parameters( filename )

    I = load( filename );
    
    [P.Ny,P.Nx] = size(I.A);
    
    % CHOOSE MODEL TYPE: 
    % TRUE  = (STERN LAYER + DIFFUSE LAYER) 
    % FALSE = (DIFFUSE LAYER)
    P.STERN = false;
    
    if ( P.STERN == true )
        P.MODEL = 'STERN';
    else
        P.MODEL = 'DIFFUSE';
    end
    
    % BOUNDARY CONDITIONS: POTENTIAL
    % CHOICE 1 or 2 or 3
    P.BC_U = 3; % 1 - Zero Dirichlet , 2,3- Neumann, non-zero Drichlet
    % IN rocks there is no change in surf. charghe density (^^)
    
    % BOUNDARY CONDITION: CATION & ANION
    % CHOSE 1 --- or 2 for zero flux bnd cdn
    P.BC_CONC = 2;  
  
    % CONSTANTS
    P.Faraday   = 95484.56;
    P.epsilonO  = 8.85*1e-12;
    P.epsilonR  = 80; % permitivity of electrolyte soln
    P.kB        = 1.38*1e-23;
    P.T         = 289.85;
    P.NA        = 6.02*1e+23;
    P.z         = 1;
    P.e         = 1.60*1e-19;
  
    P.c         = P.Faraday/(P.epsilonO*P.epsilonR); % abs. relative * rel. epsilon
   
    % SPACE DISCRETIZATION
    P.Lx        = 16*(1e-9);              % DOMAIN LENGTH
    P.Ly        = 16*(1e-9);              % DOMAIN LENGTH
    P.hx        = P.Lx*(1/(P.Nx+1));
    P.hy        = P.Ly*(1/(P.Ny+1));
  
    % TIME DISCRETIZATION
    P.Nt        = 512*512; % TIME POINTS (Stiff Problem)
    P.t_start   = 0;
    P.t_end     = 4.1e-7; % crit
    P.dt        = (P.t_end - P.t_start)/P.Nt;  
  
    P.Co        = 1e+1;                   % CONCENTRATION 1 mol/m3 or 10 mol/m3
    P.sigma     = 1e-3*0;                   % SURFACE CHARGE AT ROCK SURFACE
    P.U_ext     = 1e-3;                 % 1 V/m SCALED TO NANOMETER
       
    P.m_p       = 5.00*1e-8;              % CATION MOBILITY  COEFFICIENT        
    P.m_n       = 5.00*1e-8;              % ANION  MOBILITY  COEFFICIENT       
    P.D_p       = (P.m_p*P.kB*P.T)/(P.e); % CATION DIFFUSION COEFFICIENT         
    P.D_n       = (P.m_n*P.kB*P.T)/(P.e); % ANION  DIFFUSION COEFFICIENT
      
    % Solvers and Algorithms Types
    % ALGO 1: Simple  GS Solver (U) + Explicit Euler (Cp Cn)
    % ALGO 2: Simple SOR Solver (U) + Explicit Euler (Cp Cn)
    % ALGO 3: Complex GS Solver (U) + Explicit Euler (Cp Cn)
    % ALGO 4: Simple  GS Solver (U) + Implicit Euler (Cp Cn)
    % ALGO 5: Simple  CG Solver (U) + Explicit Euler (Cp Cn)
  
    P.ALGO      = 2;
    P.max_iter  = 1e+2;       
    P.threshold = 1e-3;   
    P.SOR_w     = 1.7;
  
    % SCALING VALUES
    P.SratioD = 10;     % STERN-DIFFUSE LENGTH RATION    
    P.Debye   = sqrt((P.epsilonO*P.epsilonR*P.kB*P.T)/(2*P.NA*((P.z*P.e)^2)*P.Co));
    P.Ut      = (P.kB*P.T)/(P.z*P.e); % thermal voltage
    P.D       = (P.D_p + P.D_n)/2; % avg. cation anion
    P.V       = P.sigma; % surf. charge density
    
    % Electrode particle properties
    P.D_s  =  P.D_p* 10^(-3);
    P.C_smax  =  P.Co * 10;
    
    % DIFFUSE LAYER  
    D.Debye   = P.Debye; 
    D.Ut      = P.Ut;
    D.D       = P.D;
      
    D.c       = P.c/(D.Debye*D.Ut*P.NA);      
     
    D.Lx      = P.Lx/D.Debye;
    D.hx      = D.Lx*(1/(P.Nx+1));
    D.Ly      = P.Ly/D.Debye;
    D.hy      = D.Ly*(1/(P.Ny+1));      
      
    D.Nt      = P.Nt;
    D.t_start = P.t_start;
    D.t_end   = (P.t_end*D.D)/(D.Debye^2); 
    D.dt      = (D.t_end - D.t_start)/P.Nt; 
      
    D.Co      = P.Co*(P.NA*(D.Debye)^3);
    D.sigma   = P.sigma/D.Ut; 
    D.U_ext   = P.U_ext/D.Ut;
        
    D.m_p     = (P.m_p*D.Ut)/D.D;        
    D.m_n     = (P.m_n*D.Ut)/D.D;       
    D.D_p     = (P.D_p)/D.D;       
    D.D_n     = (P.D_n)/D.D;
      
    % STERN LAYER
    S.Debye   = D.Debye/P.SratioD; 
    S.Ut      = D.Ut;
    S.D       = D.D/P.SratioD; % LOWER DIFFUSION & MOBILITY IN STERN LAYER

    S.c       = P.c/(S.Debye*S.Ut*P.NA);       
    
    P.m_p_stern  = (5.00*1e-8)/P.SratioD;        % CATION MOBILITY  COEFFICIENT        
    P.m_n_stern  = (5.00*1e-8)/P.SratioD;        % ANION  MOBILITY  COEFFICIENT       
    P.D_p_stern  = (P.m_p_stern*P.kB*P.T)/(P.e); % CATION DIFFUSION COEFFICIENT         
    P.D_n_stern  = (P.m_n_stern*P.kB*P.T)/(P.e); % ANION  DIFFUSION COEFFICIENT  
             
    S.Lx      = P.Lx/S.Debye;      
    S.hx      = S.Lx*(1/(P.Nx+1));  
    S.Ly      = P.Ly/S.Debye;       
    S.hy      = S.Ly*(1/(P.Ny+1));      
      
    S.Nt      = P.Nt;
    S.t_start = P.t_start;
    S.t_end   = (P.t_end*S.D)/(S.Debye^2);
    S.dt      = (S.t_end - S.t_start)/P.Nt; 
      
    S.Co      = P.Co*(P.NA*(S.Debye)^3);
    S.sigma   = P.sigma/S.Ut;
    S.U_ext   = P.U_ext/S.Ut;
        
    S.m_p     = (P.m_p_stern*S.Ut)/S.D;        
    S.m_n     = (P.m_n_stern*S.Ut)/S.D;       
    S.D_p     = (P.D_p_stern)/S.D;       
    S.D_n     = (P.D_n_stern)/S.D;   
    
end 
