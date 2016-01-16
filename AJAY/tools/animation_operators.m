
D = size(Ut);

figure(12)
hold on

for(t=1:D(3))
    
    clf
    
    [  P1 P2 P3 N1 N2 N3  ] = analysis_operators( P, F, t, Cpt, Cnt, Ut );
    
    pause(1)
    
end

