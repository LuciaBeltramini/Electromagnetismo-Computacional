% Rectangular device

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

% DEFINE PLANE WAVE SOURCE PARAMETERS
f0    = 60.0 * gigahertz;
theta = 30*degrees;
lam0  = c0/f0;
k0    = 2*pi/lam0;
MODE  = 'E';

% DEFINE DIFFRACTION GRATING PARAMETERS
L   = 1.5 * centimeters;
d   = 1.0 * centimeters;
er1 = 1.0;
er2 = 9.0;

% DEFINE FDFD PARAMETERS
NRES   = 100;
SPACER = lam0*[1 1];
NPML   = [20 20 20 20];
ermax  = max([er1 er2]);
nmax   = sqrt(ermax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID RESOLUTION
dx = lam0/nmax/NRES;
dy = lam0/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(L/dx);
dx = L/nx;

ny = ceil(d/dy);
dy = d/ny;

% GRID SIZE
Sx = L;
Nx = ceil(Sx/dx);
Sx = Nx*dx;

Sy = SPACER(1) + d + SPACER(2);
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
%% BUILD GEOMETRY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE TO AIR
ER2 = er1*ones(Nx2,Ny2);
UR2 = ones(Nx2,Ny2);

% Dimensiones del rectángulo
rect_width  = 0.2 * centimeters;
rect_height = 0.2 * centimeters;

% Centro del dominio (2x grid)
center_x = round(Nx2/2);
center_y = round(Ny2/2);

% Convertir dimensiones a celdas
half_w = round((rect_width/dx2)/2);
half_h = round((rect_height/dy2)/2);

% Índices del rectángulo
ix1 = center_x - half_w;
ix2 = center_x + half_w;
iy1 = center_y - half_h;
iy2 = center_y + half_h;

% Asignar er2 en el rectángulo
ER2(ix1:ix2, iy1:iy2) = er2;

% Mostrar geometría
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

% INCIDENT WAVE VECTOR
nsrc  = sqrt(UR2(1,1)*ER2(1,1));
kxinc = k0*nsrc*sin(theta);
kyinc = k0*nsrc*cos(theta);
kinc  = [ kxinc ; kyinc ];

% BUILD DERIVATIVE MATRICES
NS  = [Nx Ny];
RES = [dx dy];
BC  = [1 1];
[DEX,DEY,DHX,DHY] = yeeder2d(NS,k0*RES,BC,kinc/k0);
        
% BUILD WAVE MATRIX
if MODE == 'E'
    A = DHX/URyy*DEX + DHY/URxx*DEY + ERzz;
else
    A = DEX/ERyy*DHX + DEY/ERxx*DHY + URzz;
end

% CALCULATE SOURCE FIELD
fsrc  = exp(-1i*(kxinc*X + kyinc*Y));

% CALCULATE SCATTERED-FIELD MASKING MATRIX
ny = NPML(1) + 2;
Q  = zeros(Nx,Ny);
Q(:,1:ny) = 1;
Q  = diag(sparse(Q(:)));

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
