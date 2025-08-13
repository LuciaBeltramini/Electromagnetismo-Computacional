% Simulacion de onda cilindrica propagandose en el vacio.

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS
degrees     = pi/180;
meters      = 1;
centimeters = 1e-2 * meters;
millimeters = 1e-3 * meters;
micrometers = 1e-6 * meters;
nanometers  = 1e-9 * meters;
inches      = 2.54 * centimeters;
feet        = 12 * inches;
seconds     = 1;
hertz       = 1/seconds;
kilohertz   = 1e3 * hertz;
megahertz   = 1e6 * hertz;
gigahertz   = 1e9 * hertz;
terahertz   = 1e12 * hertz;
petahertz   = 1e15 * hertz;

% CONSTANTS
e0 = 8.85418782e-12 * 1/meters;
u0 = 1.25663706e-6 * 1/meters;
N0 = sqrt(u0/e0);
c0 = 299792458 * meters/seconds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE CILINDRCAL WAVE SOURCE PARAMETERS
f0    = 120.0 * gigahertz;
theta = 30*degrees;
lam0  = c0/f0;
k0    = 2*pi/lam0;
MODE  = 'E';

% MEDIUM PARAMETERS

er1 = 1.0;
ur1 = 1.0;

% DEFINE FDFD PARAMETERS
NRES   = 100;
real_x = 1.5 * centimeters;
real_y =  1.0 * centimeters;
SPACER = lam0*[1 1];
NPML   = [20 20 20 20];
ermax  = er1;
nmax   = sqrt(ermax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID RESOLUTION
dx = lam0/nmax/NRES;
dy = lam0/nmax/NRES;

% GRID SIZE
Sx = real_x;
Nx = ceil(Sx/dx);
Sx = Nx*dx;

Sy = SPACER(1) + real_y + SPACER(2);
Ny = NPML(1) + ceil(Sy/dy) + NPML(1);
Sy = Ny*dy;

% 2X GRID
Nx2 = 2*Nx;             dx2 = dx/2;
Ny2 = 2*Ny;             dy2 = dy/2;

% CALCULATE AXIS VECTORS
xa = [1:Nx]*dx;
ya = [1:Ny]*dy;
[Y,X] = meshgrid(ya,xa);
xa2 = [1:Nx2]*dx2;
ya2 = [1:Ny2]*dy2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DIFFRACTION GRATING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE VACUMM
ER2 = er1*ones(Nx2,Ny2);
UR2 = ur1*ones(Nx2,Ny2);
 
 
% SHOW REPRESENTED SPACE
subplot(121);
imagesc(xa2,ya2,ER2.');
axis equal tight off;
colorbar;
title('ER2');

% INCORPORATE PML
[ERxx,ERyy,ERzz,URxx,URyy,URzz] ...
                    = addupml2d(ER2,UR2,NPML);

% DIAGONALIZE MATERIAL TENSORS
ERxx = diag(sparse(ERxx(:)));
ERyy = diag(sparse(ERyy(:)));
ERzz = diag(sparse(ERzz(:)));
URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
URzz = diag(sparse(URzz(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FDFD ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BUILD DERIVATIVE MATRICES
NS  = [Nx Ny];
RES = [dx dy];
BC  = [0 0];
[DEX,DEY,DHX,DHY] = yeeder2d(NS,k0*RES,BC,k0/k0);
        
% BUILD WAVE MATRIX
if MODE == 'E'
    A = DHX/URyy*DEX + DHY/URxx*DEY + ERzz;
else
    A = DEX/ERyy*DHX + DEY/ERxx*DHY + URzz;
end

% CALCULATE SOURCE FIELD
x0 = Sx/2;                % posición fuente en x
y0 = Sy/2;                % posición fuente en y
R  = sqrt( (X - x0).^2 + (Y - y0).^2 );
R(R==0) = dx/10;        
fsrc = besselh(0, 1, k0 * R);  % Funcion de Hankel primera especie, orden 0

% CALCULATE SCATTERED-FIELD MASKING MATRIX
Q = zeros(Nx,Ny);
r_mask = lam0*2.5;  % radio zona TF
Q( R < r_mask ) = 1;
Q = diag(sparse(Q(:)));


% CALCULATE SOURCE VECTOR
b = (Q*A - A*Q)*fsrc(:);

% SOLVE FOR FIELD
f = A\b;
f = reshape(f,Nx,Ny);

% SHOW FIELD
subplot(122);
pcolor(xa,ya,real(f).');
shading interp;
axis equal tight;
set(gca,'YDir','reverse');
caxis([-1 1]);
colorbar;
title('real(E_z)');

%%
figure(2)
pcolor(xa,ya,abs(f).'); shading interp; axis equal tight;
set(gca,'YDir','reverse'); colorbar; title('|E_z|');


%% Visualización de la máscara Q
figure(3);
Q_matrix = reshape(full(diag(Q)),Nx,Ny);
imagesc(xa,ya,Q_matrix.');
axis equal tight;
set(gca,'YDir','reverse');
colormap(parula);
colorbar;
title('Máscara Q (zona de inyección)');
xlabel('x [m]');
ylabel('y [m]');
