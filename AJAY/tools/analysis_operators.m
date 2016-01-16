
function [  P1 P2 P3 N1 N2 N3  ] = analysis_operators( P, F, t, Cpt, Cnt, Ut )

  c       = P.c;
  sigma   = P.sigma;
  epsilon = P.epsilon;
  D_p     = P.D_p;
  D_n     = P.D_n;
  m_p     = P.m_p;
  m_n     = P.m_n;
  
  Uext    = P.Uext;
  f       = P.f;
  
  Nx  = P.Nx;
  Ny  = P.Ny;
  hx  = P.hx;
  hy  = P.hy;
  hxy = (2/(hx*hx)) + (2/(hy*hy));  
  dt  = P.dt;
  
  Cpo(:,:) = Cpt(:,:,t); Cno(:,:) = Cnt(:,:,t); Un(:,:) = Ut(:,:,t);
  
  P1 = zeros(P.Nx+2,P.Ny+2);
  P2 = zeros(P.Nx+2,P.Ny+2);
  P3 = zeros(P.Nx+2,P.Ny+2);  
  
  N1 = zeros(P.Nx+2,P.Ny+2);
  N2 = zeros(P.Nx+2,P.Ny+2);
  N3 = zeros(P.Nx+2,P.Ny+2);  
  
  for j = 2:Ny+1
    
    for i = 2:Nx+1
      
      if (F.FLAG(j,i) >= F.C_F)

        % CATION
        % CALCULATING TERM: P1
        P1(j,i) = D_p*(  ((Cpo(j,i+1) - 2*Cpo(j,i) + Cpo(j,i-1))/(hx*hx)) ...
                       + ((Cpo(j+1,i) - 2*Cpo(j,i) + Cpo(j-1,i))/(hy*hy)) );
                   
        % CALCULATING TERM: P2
        P2(j,i) = m_p*(  ((Cpo(j,i+1) - Cpo(j,i-1))*(Un(j,i+1) - Un(j,i-1))/(4*hx*hx)) ...
                       + ((Cpo(j+1,i) - Cpo(j-1,i))*(Un(j+1,i) - Un(j-1,i))/(4*hy*hy)) );
                    
        % CALCULATING TERM: P3
        P3(j,i) = m_p*Cpo(j,i)*(  ((Un(j,i+1) - 2*Un(j,i) + Un(j,i-1))/(hx*hx)) ...
                                + ((Un(j+1,i) - 2*Un(j,i) + Un(j-1,i))/(hy*hy)) );
       
         
        % ANION
        % CALCULATING TERM: N1
        N1(j,i) = D_n*(  ((Cno(j,i+1) - 2*Cno(j,i) + Cno(j,i-1))/(hx*hx)) ...
                       + ((Cno(j+1,i) - 2*Cno(j,i) + Cno(j-1,i))/(hy*hy)) );
                   
        % CALCULATING TERM: N2
        N2(j,i) = - m_n*(  ((Cno(j,i+1) - Cno(j,i-1))*(Un(j,i+1) - Un(j,i-1))/(4*hx*hx)) ...
                         + ((Cno(j+1,i) - Cno(j-1,i))*(Un(j+1,i) - Un(j-1,i))/(4*hy*hy)) );
                    
        % CALCULATING TERM: N3
        N3(j,i) = - m_n*Cno(j,i)*(  ((Un(j,i+1) - 2*Un(j,i) + Un(j,i-1))/(hx*hx)) ...
                                  + ((Un(j+1,i) - 2*Un(j,i) + Un(j-1,i))/(hy*hy)) );
                                                           
      end
    
    end

  end
  
    PS = P.dt*(P1 + P2 + P3);
    NS = P.dt*(N1 + N2 + N3);  
    
    Cpn = Cpo + PS;
    Cnn = Cno + NS;
  
    x = (P.Lx-1)*(0:1/(P.Nx-1):1); y = (P.Lx-1)*(0:1/(P.Ny-1):1);
    [X Y] = meshgrid(x,y);
    
    figure(12)
    hold on
    grid on
    
    subplot(2,5,1)
    surf(X,Y,P1([2:P.Nx+1],[2:P.Ny+1]))  
    title('+ D_p \Delta C_p')    
    subplot(2,5,2)
    surf(X,Y,P2([2:P.Nx+1],[2:P.Ny+1]))  
    title('+ \mu_p \nabla C_p \nabla U')    
    subplot(2,5,3)
    surf(X,Y,P3([2:P.Nx+1],[2:P.Ny+1]))  
    title('+ \mu_p C_p \Delta U')  
    subplot(2,5,4)
    surf(X,Y,PS([2:P.Nx+1],[2:P.Ny+1]))  
    title('+ \Sigma T_p')    
    subplot(2,5,5)
    surf(X,Y,Cpn([2:P.Nx+1],[2:P.Ny+1]))  
    title('C_p')     
    
    subplot(2,5,6)
    surf(X,Y,N1([2:P.Nx+1],[2:P.Ny+1]))  
    title('+ D_n \Delta C_n')   
    subplot(2,5,7)
    surf(X,Y,N2([2:P.Nx+1],[2:P.Ny+1]))  
    title('- \mu_n \nabla C_n \nabla U')    
    subplot(2,5,8)
    surf(X,Y,N3([2:P.Nx+1],[2:P.Ny+1]))       
    title('- \mu_n C_n \Delta U')
    subplot(2,5,9)
    surf(X,Y,NS([2:P.Nx+1],[2:P.Ny+1]))  
    title('+ \Sigma T_n')   
    subplot(2,5,10)
    surf(X,Y,Cnn([2:P.Nx+1],[2:P.Ny+1]))  
    title('C_n')     
    
end
