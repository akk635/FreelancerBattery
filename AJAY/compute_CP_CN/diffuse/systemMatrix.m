function [ A ] = systemMatrix(  P, S, D, F, Un )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    A = sparse(P.Nx+2 * P.Ny+2, P.Nx+2 * P.Ny+2);
    
    c     = D.c;
    sigma = D.sigma;
    D_p   = D.D_p;
    D_n   = D.D_n;
    m_p   = D.m_p;
    m_n   = D.m_n;
  
    Nx  = P.Nx;
    Ny  = P.Ny;
    hx  = D.hx;
    hy  = D.hy;
    hxy = (2/(hx*hx)) + (2/(hy*hy));  
    dt  = D.dt;
    DNx = P.Nx+2;
    DNy = P.Ny+2;

    for index = 1 : DNx * DNy
        rowIdx = index / DNx; % Domain indices
        colIdx = mod(index , DNx);
        % Boundary Conditions
        if (rowIdx == 0) % north boundary
            A(index, index) = -1;
            A(index, ((rowIdx + 1) * DNx + colIdx) = 1;
        elseif (rowIdx == DNy - 1) % south boundary
            A(index, index) = 1;
            A(index , ((rowIdx - 1) * DNx + colIdx) = -1;
        elseif (colIdx == 1) % west boundary
            A(index, index) = -1;
            A(index, index + 1) = 1;
        elseif (colIdx == 0)
            A(index, index) = 1;
            A()
                    
            
        A (index, (rowIdx * P.Nx + 2) + colIdx ) = 0;     
    end
     
            
end

