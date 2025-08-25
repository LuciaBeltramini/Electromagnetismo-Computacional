% 1- Codigo para calcular la tramitancia y la refelctancia de una onda electromagnetica al cambiar de medio, siendo este SILICIO NANOPOROSO. 

% 2- Se necesita LA FUNCION 'indiceSi5.m' que contengan el indice de refraccion del Si cristalino para distintas longitudes de onda.

% 3- Este codigo solo funciona para incidencia normal y calculando E mode.

% 4 - Se da la opcion de colocar rugosidades.


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
lam2 = 800 * nanometers;
delta = 10 * nanometers;
wavelength = lam1:delta:lam2;
rep = 80;

% FOR SAVING IMPORTANT DATA
T = [];
R = [];
AA = [];

% SWEEP
for lam0 = wavelength
    fprintf('------------------ LAMDA [nm] = %.1f ------------------\n',  lam0/nanometers); 
    t_rep = [];
    r_rep = [];
    a_rep = [];

    for times=1:rep
        fprintf('Repeticion n° = %.1f\n',  times); 
        % REFRACTIVE INDEX OF Si AND AIR
        n_Si = indiceSi5(lam0);
        
        % AIR
        n_air = 1;
        
        % EFFECTIVE REFRACTIVE INDEX
        p = 0.46; 
        n_eff = ((n_Si.^(2/3)) * (1 - p) + (n_air.^(2/3)) * p).^(3/2);
    
        % DEFINE PLANE WAVE SOURCE PARAMETERS
        theta = 0*degrees;
        f0 = c0/lam0;
        k0    = 2*pi/lam0;
        MODE  = 'E';
        
        % DEFINE MEDIA PARAMETERS
        er1 = 1.0;
        er2 = n_eff^2;
        width  = 2.3960e-06;
        height = 2.1994e-06;
    
        % ADD ROUGHNESS TO THE INTERFACE
        roughness = true;
        sigma = 4e-8;  
    
        % DEFINE FDFD PARAMETERS
        NRES   = 25;
        SPACER = lam0*[4 4];
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
        
        Sy = SPACER(1) + height + SPACER(2);
        Ny = NPML(1) + ceil(Sy/dy) + NPML(2);
        Sy = Ny*dy;
        
        % ZONE OF INTEGRATION (ZI) FOR R & T
        y_ref = round(Ny/2) - round(SPACER(1)/1.5/dy);
        y_trn = round(Ny/2) + round(SPACER(2)/1.5/dy);
        
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
    
        if roughness == true
            height = height + randn*sigma;         
        end

        h = round(height/2/dy2);
        ER2(:, center_y - h:center_y + h) = er2;
    
        % SHOW BOTH MEDIUMS AND ZI
        if length(wavelength)==1 & times ==1
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
            colormap(flipud(cool));
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
        if length(wavelength)==1 & times==1
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
            colormap(flipud(cool));
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
            Sy_inc_line  = 0.5*real( Ez_inc(:, y_ref)  .* conj(Hx_inc(:,  y_ref)/urref) );  
            Sy_ref_line  = 0.5*real( Ez_scat(:,y_ref)  .* conj(Hx_scat(:, y_ref)/urref) );
            Sy_trn_line  = 0.5*real( Ez_tot(:, y_trn)  .* conj(Hx_tot(:,  y_trn)/urtrn) );
            
            % POWER
            P_inc = trapz(xa, Sy_inc_line);          
            P_ref = -trapz(xa, Sy_ref_line);         
            P_trn =  trapz(xa, Sy_trn_line);  
            
            % TRANSMITTANCE AND REFRECTANCE
            r = real(P_ref / P_inc);
            t = real(P_trn / P_inc);
            a = 1 - t - r;
            
            r_rep(end+1) = r;
            t_rep(end+1) = t;
            a_rep(end+1) = a;
        
        end 
    end

    R(end+1) = mean(r_rep);
    T(end+1) =  mean(t_rep);
    AA(end+1) = mean(a_rep);

    fprintf('Reflectancia  R = %.4f\n',  mean(r_rep));
    fprintf('Transmitancia T = %.4f\n',  mean(t_rep));
    fprintf('Absorbancia A = %.4f\n',  mean(a_rep)); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT TRANSMITTANCE AND REFRECTANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(wavelength) ~= 1
    figure(10)
    hold on
    plot(wavelength * 1/nanometers, R*100, 'o-', 'LineWidth', 1); hold on;
    plot(wavelength * 1/nanometers, T*100, 'o-', 'LineWidth', 1); hold on;
    plot(wavelength * 1/nanometers, AA*100,'o-', 'LineWidth', 1)
    xlabel('Longitud de onda [nm]');      
    ylabel('R y T porcentual'); 
    title('Reflectancia y Transmitancia vs Longitud de onda');
    grid on;
    ylim([0 105]);
    legend({'R (Reflectancia)','T (Transmitancia)', 'A (Absrobancia)'}, 'Location','best');   
end