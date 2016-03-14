%% MESSAGE TO CHOOSE DIMENSION OR DIMENSIONLESS
scale = 1; 
if ( scale == 1)
    rescaleConc = (P.NA*(D.Debye)^3);
    rescaleU    = P.Ut; 
    L           = P.Lx;
else
    rescaleConc = 1;
    rescaleU    = 1;
    L           = 1;
end

%% SCALING CATION AND ANION CONCENTRATIONS
cation = Cpt([2:end-1],[2:end-1],:)/rescaleConc;
anion  = Cnt([2:end-1],[2:end-1],:)/rescaleConc;
en     = cation - anion;
% FINDING MIN & MAX CATION & ANION CONCENTRATION
MIN_CO = min( min(min(min(cation(cation>0)))), min(min(min(anion(anion>0)))) );
MAX_CO = max( max(max(max(cation))) , max(max(max(anion))) );      
MIN_EN = min(min(min(en)));  
MAX_EN = max(max(max(en)));   

%% SCALING POTENTIAL
potential = Ut([2:end-1],[2:end-1],:)*rescaleU;
% FINDING MIN & MAX POTENTIAL
MIN_U = min(min(min(potential)));  
MAX_U = max(max(max(potential)));

%% THE SPATIAL GRID
x = linspace(0,1*L,size(potential,1)+1);
y = linspace(0,1*L,size(potential,1)+1);
[XX,YY] = meshgrid(x,y);
ZZ = zeros(size(XX));

%% MASKING ALL THE OBJECT CELLS
XC = XX(1:end-1,1:end-1);
YC = YY(1:end-1,1:end-1);
dx = XX(1,2)-XX(1,1);
dy = YY(2,1)-YY(1,1);
XC = XC + dx/2;
YC = YC + dy/2;

% MASKING ANION
CMASK_A = anion(:,:,end);
CMASK_A(~F.GEO([2:end-1],[2:end-1])) = NaN;

% MASKING CATION
CMASK_C = cation(:,:,end);
CMASK_C(~F.GEO([2:end-1],[2:end-1])) = NaN;

% MASKING ION DIFFERENCE
CMASK_EN = en(:,:,end);
CMASK_EN(~F.GEO([2:end-1],[2:end-1])) = NaN;

% MASKING POTENTIAL
CMASK_U = potential(:,:,end);
CMASK_U(~F.GEO([2:end-1],[2:end-1])) = NaN;

figure(1);
subplot(2,2,1)
levels_U = linspace(MIN_U,MAX_U,(MAX_U-MIN_U)+1);
% THE OFFSET IS NEEDED FOR PLOTTING PURPOSES
offset = abs(floor( min(min(CMASK_U))));
surf(XX,YY,ZZ-offset,CMASK_U);
hold on
contour3(XC,YC,CMASK_U,'k');
%[Co,ch] = contour3(XC,YC,CMASK_U,levels_U,'k');
%clabel(Co,ch,'FontSize',9,'FontWeight','demi');
shading flat
axis equal
view(2)
box on
grid off
colorbar
set(gca,'Layer','top');
title('POTENTIAL [V]')
xlabel('Distance [m]')
ylabel('Distance [m]')

%figure(2)
subplot(2,2,2)
levels_EN = linspace(MIN_EN,MAX_EN,(MAX_EN-MIN_EN)+1);
% THE OFFSET IS NEEDED FOR PLOTTING PURPOSES
offset = abs(floor(min(min(CMASK_EN))));
surf(XX,YY,ZZ-offset,CMASK_EN);
hold on
contour3(XC,YC,CMASK_EN,'k');
shading flat
axis equal
view(2)
box on
grid off
colorbar
set(gca,'Layer','top');
title('EN = CATION - ANION  [mol/m^3]')
xlabel('Distance [m]')
ylabel('Distance [m]')

%figure(3)
subplot(2,2,3)
levels_CO = linspace(MIN_CO,MAX_CO,(MAX_CO-MIN_CO)+1);
% THE OFFSET IS NEEDED FOR PLOTTING PURPOSES
offset = abs(floor(min(min(CMASK_C))));
%[Co,ch] = contour3(XC,YC,CMASK_C,levels_CO,'k');
%clabel(Co,ch,'FontSize',9,'FontWeight','demi');
surf(XX,YY,ZZ-offset,CMASK_C);
hold on
contour3(XC,YC,CMASK_C,'k');
shading flat
axis equal
view(2)
box on
grid off
colorbar
set(gca,'Layer','top');
title('CATION CONC. [mol/m^3]')
xlabel('Distance [m]')
ylabel('Distance [m]')

%figure(4)
subplot(2,2,4)
levels_CO = linspace(MIN_CO,MAX_CO,(MAX_CO-MIN_CO)+1);
% THE OFFSET IS NEEDED FOR PLOTTING PURPOSES
offset = abs(floor(min(min(CMASK_A))));
%[Co,ch] = contour3(XC,YC,CMASK_A,levels_CO,'k');
%clabel(Co,ch,'FontSize',9,'FontWeight','demi');
surf(XX,YY,ZZ-offset,CMASK_A);
hold on
contour3(XC,YC,CMASK_A,'k');
shading flat
axis equal
view(2)
box on
grid off
colorbar
set(gca,'Layer','top');
title('ANION CONC. [mol/m^3]')
xlabel('Distance [m]')
ylabel('Distance [m]')
