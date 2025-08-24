% Codigo para calcular la tramitancia y la refelctancia de una onda
% electromagnetica al cambiar de medio. 
% Este codigo solo funciona para incidencia normal y calculando E mode.


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

% WAVELENGTH (EMPTY SPACE)
lam1 = 700 * nanometers;
lam2 = 700 * nanometers;
delta = 100 * nanometers;
wavelength = lam1:delta:lam2;

% FOR SAVING IMPORTANT DATA
T = [];
R = [];

% DEFINE PLANE WAVE SOURCE PARAMETERS
for lam0 = wavelength
    theta = 0*degrees;
    f0 = c0/lam0;
    k0    = 2*pi/lam0;
    MODE  = 'E';
    
    % DEFINE DIFFRACTION GRATING PARAMETERS
    er1 = 1.0;
    er2 = 5.0;
    
    % DIMENSIONS OF THE GRID
    width  = lam0*6;
    height = lam0*6;
    
    % DEFINE FDFD PARAMETERS
    NRES   = 40;
    NPML   = [20 20];
    ermax  = max([real(er1) real(er2)]);
    nmax   = sqrt(ermax);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CALCULATE OPTIMIZED GRID
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % GRID RESOLUTION
    dx = lam0/nmax/NRES;
    dy = lam0/nmax/NRES;
    
    % SNAP GRID TO CRITICAL DIMENSIONS
    nx = ceil(width/dx);
    dx = width/nx;
    
    ny = ceil(height/dy);
    dy = height/ny;
    
    % GRID SIZE
    Sx = width;
    Nx = ceil(Sx/dx);
    Sx = Nx*dx;
    
    Sy = height;
    Ny = NPML(1) + ceil(Sy/dy) + NPML(2);
    Sy = Ny*dy;
    
    % ZONE OF INTEGRATION (ZI) FOR R & T
    y_ref = round(Ny/2) - round(2*lam0/dy);
    y_trn = round(Ny/2) + round(2*lam0/dy);
    
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
    
    
    % BUILD RECTANGLE INTO X2 GRID
    center_y = round(Ny2/2);
    
    % CHANGE ER
    ER2(:, center_y:end) = er2;
    
    % SHOW BOTH MEDIUMS AND ZI
    if length(wavelength)==1
        subplot(121);
        imagesc(xa2,ya2,real(ER2).');
        axis equal tight;
        colorbar;
        xlabel('Eje X [m]');
        ylabel('Eje Y [m]');
        title('Er DEL MEDIO');
        hold on;
        y_ref_pos = ya(y_ref);   
        y_trn_pos = ya(y_trn);   
        line([min(xa2) max(xa2)], [y_ref_pos y_ref_pos], 'Color','k','LineWidth',1.5);
        line([min(xa2) max(xa2)], [y_trn_pos y_trn_pos], 'Color','k','LineWidth',1.5);
    hold off;
    end
    
    % INCORPORATE PML
    [ERxx,ERyy,ERzz,URxx,URyy,URzz] ...
                        = addupml2d(ER2,UR2,[0 0 NPML]);
    
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
    BC  = [1 0];
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
    if length(wavelength)==1
        subplot(122);
        pcolor(xa,ya,real(f).');
        shading interp;
        axis equal tight;
        set(gca,'YDir','reverse');
        caxis([-1 1]);
        colorbar;
        xlabel('Eje X [m]');
        ylabel('Eje Y [m]');
        title('R(Ez)');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TRANSMITTANCE AND REFRECTANCE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if theta==0
        % EXTRACT MATERIAL PROPERTIES IN EXTERNAL REGIONS
        urref = UR2(1,1);
        urtrn = UR2(1,Ny2);
        erref = ER2(1,1);
        ertrn = ER2(1,Ny2);
        nref  = sqrt(urref*erref);
        ntrn  = sqrt(urtrn*ertrn);
        
        % TOTAL, SOURCE AND SCATTERED FIELDS
        Ez_tot  = f;        
        Ez_inc  = fsrc;     
        Ez_scat = Ez_tot - Ez_inc;
        
        % PARTIAL DERIVATIVES OF Ez
        [Ezy_tot, Ezx_tot]   = gradient(Ez_tot, dy, dx); 
        [Ezy_inc, Ezx_inc]   = gradient(Ez_inc, dy, dx);
        [Ezy_scat,Ezx_scat]  = gradient(Ez_scat,dy, dx);
        
        %CALCULATE MAGNETIC FIELDS
        omega = 2*pi*f0;
        Hx_tot  = -(1./(1i*omega*u0)) .* Ezy_tot;
        Hx_inc  = -(1./(1i*omega*u0)) .* Ezy_inc;
        Hx_scat = -(1./(1i*omega*u0)) .* Ezy_scat;
        
        % POYTING FLUX IN EACH ZI
        Sy_inc_line  = 0.5*real( Ez_inc(:, y_ref)  .* conj(Hx_inc(:,  y_ref)/urref) );    %agrego los ur_xx de cada medio
        Sy_ref_line  = 0.5*real( Ez_scat(:,y_ref)  .* conj(Hx_scat(:, y_ref)/urref) );
        Sy_trn_line  = 0.5*real( Ez_tot(:, y_trn)  .* conj(Hx_tot(:,  y_trn)/urtrn) );
        
        % POWER
        P_inc = trapz(xa, Sy_inc_line);          
        P_ref = -trapz(xa, Sy_ref_line);         
        P_trn =  trapz(xa, Sy_trn_line);  
        
        % TRANSMITTANCE AND REFRECTANCE
        r = real(P_ref / P_inc);
        t = real(P_trn / P_inc);
        
        R(end+1) = r;
        T(end+1) = t;
    
        fprintf('Reflectancia  R = %.4f\n', r);
        fprintf('Transmitancia T = %.4f\n', t);
        fprintf('Verifico R+T = %.4f\n', r+t); 
    
    end 
end
% PLOT TRANSMITTANCE AND REFRECTANCE
if length(wavelength) ~= 1
    figure(1)
    plot(wavelength, R, 'o-', 'LineWidth', 2); hold on;
    plot(wavelength, T, 'o-', 'LineWidth', 2);
    xlabel('Longitud de onda [nm]');      
    ylabel('Coeficientes de Potencia'); 
    title('Reflectancia y Transmittancia vs Longitud de onda');
    grid on;
    ylim([0 1.5]);
    legend({'R (Reflectancia)','T (Transmittancia)'}, 'Location','best');
end

