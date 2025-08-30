function []=Multicapas_C()
%Diseno de un multicapa de materiales arbitrarios. Estos materiales deben
%poseer una funcion que devuelva el indice de refraccion (complejo) en
%funcion de la longitud de onda. 
%En esta versión se genera un solo archivo con todas las funciones
%Se piden una serie de datos iniciales que configuran el tipo de calculo a
%realizar, que puede ser uno o varios de los siguientes:
%1_ Calcular la reflectancia, transmitancia y absorbancia de la multicapa
%2_ Calcular el campo electrico dentro de las multicapas.
%3_ Dibujar la multicapa en diagramas de colores.
%Si se da una lista de longitudes de onda y un solo valor de ángulo
%calcula el espectro para ese valor de ángulo.
%Si se dan varios valores de angulo y una sola longitud de onda se
%dibuja en función de los primeros.
%Si se dan dos vectores, uno para el ángulo y otro para la longitud de onda construye
%gráficos de contornos. En este caso no se dibuja el campo electrico.
%Esta version utiliza la función "Rugosear3", esta permite introducir un 
%parámetro "delta" de rugosidad las interfases de la multicapa, puede ser la 
%misma para todas las interfases o tener diferentes valores. Esta Rugosidad
%se simula realizando un cambio suave de indice de refraccion en la
%interfase.
%En esta version se incorpora la posibilidad de agregar una rugosidad
%"corta" que corresponde a analizar posibles variaciones en los espesores
%de cada capa. Cada capa presenta un "delta_c" que indica la desviacion
%estandar de la fluctuación. Puede ser la misma para todas las interfases o
%tener distintos valores. 
%falta incorporar este calculo para angulos distintos.
%version 6/7/2023 Cuchu 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETROS DE INTERES PARA COMPARAR ESTE CÓDIGO CON EL DE FDFD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Polarización de la onda incidente:
p=0;               %Si es una onda s => p=0, si es una onda p => p=1, valores intermedios indican polarizaciones intermedias

%Angulo de incidencia:
tetai=0;

%Rango de valores de longitud de onda en metros:
lambdai=700e-9;
lambdaf=800e-9;
Npl=11;            %Numero de valores de longitud de onda tomados en este rango (si Npl=1 solo se toma lambdai).

%Medios:
Medios=[1 101 1];

%Tamaño de cada medio:
d=-ones(1,length(Medios)-2);
d(1)=2.1994e-06;

%Parametros de la rugosidad:
rugosidad_c=1;
delta_c=4e-8;
Nprom=200;

%Porosidad de la capa de Si:
P1=0.46;       %Si rugosidad_c es igual a 1 sortea varias configuraciones posibles a los espesores de las capas y calcula un valor medio de la relfectancia medida.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CUCHU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Definicion del problema a resolver  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%p=0;               %polarización de la onda incidente (si es una onda s => p=0, si es una onda p => p=1, valores intermedios indican polarizaciones intermedias)
lambda0=500*0.9702e-9;     %longitud de onda de referencia (en metros) que se utiliza para calcular espesores de capas como espesores ópticos a esa longitud de onda y para realizar el dibujo de la multicapa (el índice de refracción de cada capa se calcula a esa longitud de onda)

%rango de valores de longitud de onda en metros donde se realiza el cálculo del espectro
% lambdai=700e-9;
% lambdaf=800e-9;
% Npl=11;            %Numero de valores de longitud de onda tomados en este rango (si Npl=1 solo se toma lambdai).

%rango de valores de angulo en grados donde se desea calcular
%tetai=0;
tetaf=70;
Npt=1;              %Numero de valores tomados en este rango (Si Npt=1 solo realiza el cálculo para tetai).


CalculodeCampo=0;   %Si se elige igual a 1 se calcula y grafica el campo en el interior (funciona sólo si Npt=1 (un solo ángulo y rugosidad=0 y rugosidad_c=0)
ndiv=20;            %Número de divisiones en cada capa para el calculo del campo electrico dentro de la multicapa (solo se aplica si CalculodeCampo=1).


GraficaMulticapa=0; %Se elige igual a 1 si se desea que se grafique la multicapa que se pretende calcular


%Ahora debemos construir el vector de medios: Por cada medio diferente se
%debe ingresar un número diferente asociado a una funcion que devuelve el
%indice de refraccion (complejo) como función de la longitud de onda (ver mas adelante
%las definiciones precargadas con las que se cuenta)

%Medios=[1 7 101 102 101 102 101 102  101 102 101 102 101 102 101 102  101 102 101 102 101 102 101 102 101 102 101 102 101 102 101 102 101 102 101 4  ];
%Medios=[1 101 1];

%Ahora debemos ingresar el vector de espesores "d" (en metros) que tiene dos elementos menos
%que el vector de índices de refraccion (o que el vector de medios) ya que
%el primer medio corresponde al medio incidente (medio semiinfinito a la
%izquierda digamos) y el ultimo medio corresponde al sustrato (medio
%semiinfinito sobre el que está construida la multicapapa, a la derecha
%digamos). %Si el valor de espesor se coloca con signo negativo, se utiliza
%allí un valor tal que produzca un espesor óptico de 1/4 de 
%la longitud de onda de referencia (el programa calcula automaticamente
%cual debe ser el espesor de la capa correspondiente). 

% d=-ones(1,length(Medios)-2);
% 
% d(1)=2.1994e-06;
%d(1)=33e-9;
%d(2)=66e-9-25.5e-9;

% d=[3e-6 5e-6];        
%d=load('d.dat'); 

% alfa=0.99936;
% vec1 = 1:24;
% % d=[69.54 108.39]*1e-9;

% d = [ d(1) d(2) d(1) d(2) d(1) d(2) d(1) d(2) d(1) d(2) d(1) d(2) d(2) d(1) d(2) d(1) d(2) d(1) d(2) d(1) d(2) d(1) d(2) d(1) d(2)].*([1 alfa.^vec1]);


Eliminacapas=0;       %Opcion que une capas iguales contiguas para disminuir el tiempo de calculo y la reemplaza por una sola capa del espesor de todas (por default en 1).

rugosidad=0;          %Si rugosidad es igual a 1 aplica el suavizado a las interfases que corresponde a intercalar en cada interfase una serie de capas que "suavizan" el salto de indice de refraccion en la interfase
%rugosidad_c=1;        %Si rugosidad_c es igual a 1 sortea varias configuraciones posibles a los espesores de las capas y calcula un valor medio de la relfectancia medida.
 
delta=6e-8;        %Rugosidades para cada interfase (en metros). Corresponde al tamaño de la zona en que se produce el suavizado de un indice al siguiente. Tiene un elemento más que el vector de espesores, o sólo tiene un valor en cuyo caso todos los deltas se toman iguales. Este valor no se utiliza si rugosidad=0; 
%delta_c=4e-8;   %Rugosidades para cada interfase (en metros). Corresponde a la variacion del espesor de cada capa (desciacion estandar). Tiene la misma cantidad de elementos que el vector de espesores, o sólo tiene un valor en cuyo caso todos los deltas se toman iguales. Este valor no se utiliza si rugosidad_c=0; 

Nepstein=10; %número de capas usadas para el suavizado. Este valor no se utiliza si rugosidad=0; 
Ndelta=2;     %numero de veces delta que se discretiza a ambos lados de la interfase. Este valor no se utiliza si rugosidad=0; 

%Nprom=200;    %numero de evaluaciones para calcular la respuesta promedio de las variaciones de espesor en cada capa. Este valor no se utiliza si rugosidad_c=0; 

%Ahora debemos asignarle a cada medio una funcion indice de refracción

%Se pueden definir aquí valores de porosidad que se utilicen en las capas
%internas de la multicapa
%P1=0.46;
P2=0.82;
P3=1;

%o bien valores de indice que puedan ser introducidos en forma sencilla
%utilizando la funcion dummy. (ver definicion mas adelante)
n1=0.1;
n2=2;
n3=1.5;

%en lo que sigue se muestran las funciones de indice de refraccion
%incluidas en este programa. Si no se encuentra la necesaria se puede
%agregar.

Nfun(1)=cellstr('aire(lambda)');            %indice de refraccion del aire
Nfun(2)=cellstr('indiceSi5(lambda)');       %indice de refraccion del Silicio
Nfun(3)=cellstr('aluminio(lambda)');        %indice de refraccion del aluminio
Nfun(4)=cellstr('BK7(lambda)');             %indice de refraccion del vidrio mas comunmente usado en optica
Nfun(5)=cellstr('dummy(lambda,n1)');        %indice constante igual a n1
Nfun(6)=cellstr('dummy(lambda,n2)');        %indice constante igual a n2

%Nfun(7)=cellstr('dummy(lambda,n3)');        %indice constante igual a n2

Nfun(7)=cellstr('Ag(lambda)');              %indice de refraccion de la plata
Nfun(8)=cellstr('Au(lambda)');              %indice de refraccion del Oro
Nfun(9)=cellstr('Pt(lambda)');              %indice de refraccion del Platino
Nfun(10)=cellstr('PZT(lambda)');            %indice de refraccion del PZT con que se construyen tipicamente los piezoelectricos
Nfun(11)=cellstr('PMMA(lambda)');           %indice de refraccion del Polimetil-metacrilato (Acrilico)
Nfun(12)=cellstr('SnO2(lambda)');           %indice de refraccion del Oxido de Estaño
Nfun(13)=cellstr('Al2O3(lambda)');          %indice de refraccion de la alumina
Nfun(14)=cellstr('ZnO(lambda)');            %indice de refraccion del oxido de Zinc
Nfun(15)=cellstr('Siamorfo(lambda)');       %indice de refraccion del silicio amorfo
Nfun(16)=cellstr('indicePET(lambda)');      %indice de refraccion del Polietilentereftalato
Nfun(17)=cellstr('indiceSiO2(lambda)');     %indice de refraccion del Oxido de silicio (Cuarzo))
Nfun(18)=cellstr('indiceaguaLyT(lambda,T)');%indice de refraccion del agua tambien como funcion de la temperatura
Nfun(19)=cellstr('indiceagua(lambda)');     %indice de refraccion del agua
Nfun(20)=cellstr('indicealcohol(lambda)');  %indice de refraccion del alchool etilico
Nfun(21)=cellstr('indiceAIP(lambda)');      %indice de refraccion del alcohol Isopropilico
Nfun(22)=cellstr('Al2O3amorfilm(lambda)');  %indice de refraccion de la alumina amorfa depositada por sputtering
Nfun(23)=cellstr('TiO2film(lambda)');       %indice de refraccion del Oxido de titanio (depositado en forma de film)
Nfun(24)=cellstr('MAPI_Forouhi(lambda)');   %indice de refraccion de una perovkita paper con Caram
Nfun(25)=cellstr('PbI2(lambda)');           %indice de refraccion del Ioduro de Plomo
Nfun(26)=cellstr('EVA(lambda)');            %indice de refraccion del Etil Vinil Acetato
Nfun(27)=cellstr('EVA_UV(lambda)');         %indice de refraccion del Etil Vinil Acetato con proteccion UV para aumentar transmitancia en paneles solares

%tambien se pueden utilizar modelos de medio efectivo o reglas de mezcla
%para construir materiales compuestos o poroso. Debajo hay algunos ejemplos
%de las funciones que ya estan incorporadas en este mismo programa.

Nfun(101)=cellstr('looyenga(indiceSi5(lambda),1,P1)');            %Silicio Poroso (mezcla con aire) de porosidad P1. Usa el modelo de Looyenga-Lifshitz-Landau
Nfun(102)=cellstr('looyenga(indiceSi5(lambda),1,P2)');            %Silicio Poroso (mezcla con aire) de porosidad P2. Usa el modelo de Looyenga-Lifshitz-Landau
Nfun(103)=cellstr('looyenga(indiceSi5(lambda),indicealcohol(lambda),P1)');   %Silicio Poroso mojado con alcohol de porosidad P1. Usa el modelo de Looyenga-Lifshitz-Landau
Nfun(104)=cellstr('looyenga(indiceSi5(lambda),indiceagua(lambda),P1)');      %Silicio Poroso mojado con agua de porosidad P1. Usa el modelo de Looyenga-Lifshitz-Landau
Nfun(105)=cellstr('mezclahomogenea(f1,n1,n2)');                  %indice de refracción de una mezcla homogenea usando el modelo de Lorentz-Lorentz
Nfun(106)=cellstr('looyengacil(PMMA(lambda),1,P1)');             %Acrilico poroso de porosidad P1. Usa un modelo del Looyenga considerando poros cilindricos.
Nfun(107)=cellstr('Bruggeman(indiceSi5(lambda),n1,P1)');        	%Silicio Poroso lleno con un material de indice n1 con porosidad P1 . Se usa el modelo de Bruggeman.
Nfun(108)=cellstr('looyenga3m(indiceSi5(lambda),n1,n2,P1,P2)');   %Silicio Poroso lleno con dos materiales de indices n1 y n2 con fracciones en volumen P1 y P2. Se usa el modelo de Bruggeman.
Nfun(109)=cellstr('maxwell(indiceSi5(lambda),n1,P1)');        	%Silicio Poroso lleno con un material de indice n1 con porosidad P1. Se usa el modelo de Maxwell.

%hacer una funcion con el paper de Guido:Vol. 14, No. 4 / 1 Apr 2024 / Optical Materials Express
Nfun(110)=cellstr('maxwell_cil(indiceSi5(lambda),n1,P1)');        	%Silicio Poroso lleno con un material de indice n1 con porosidad P1. Se usa el modelo de Maxwell.


%%%%%%%%%%%%%%%%%%%% COMIENZA EL CÁLCULO %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
if CalculodeCampo==0
    ndiv=1;
end

Nmed=size(Medios,2);

%variable auxiliar (le agrego la longitud de onda de referencia al final
%para poder calcular los indices de refraccion en esa longitud de onda)
if Npl>1
    lambda=[lambdai:(lambdaf-lambdai)/(Npl-1):lambdaf lambda0];
else
    lambda=[lambdai lambda0];
end

%calculo el indice de refraccion para cada capa de diferente composición
aux=sort(Medios);
aux=aux([1 1+find(abs(diff(aux))>0)]);
for j=1:length(aux)
    N(aux(j),:)=eval(char(Nfun(aux(j))));
    %a=eval(char(Nfun(aux(j))));
end


%construyo el vector de espesores que tiene en cuenta si se elige un valor
%de espesor óptico de 1/4 de la longitud de onda de referencia (el indice
%de refracion a esta long de onda es calculado en el último punto del vector N)
%Este calculo lo realiza para cada espesor de capa que sea negativo (es la
%manera de indicarle al programa que ese espesor debe ser calculado como
%1/4 de lambda0.
for k=2:Nmed-1
    if d(k-1)<0
    d(k-1)=lambda0/4/real(N(Medios(k),Npl+1));
    end
end


% en el caso de que se utilice un solo valor para todas las rugosidades
% entonces se construye un vector con valores constantes. 
if length(delta)==1
    delta=delta*ones(1,length(d)+1);
end

%el programa no realiza el calculo si la longitud del vector delta no es correcta
if length(delta)~=length(Medios)-1
    disp (['Error, la cantidad de elementos en "delta" (' num2str(length(delta)) ') debe ser igual a uno o a la longitud de "d" (' num2str(length(d)) ') mas uno.'])
    return
end

%En esta seccion elimina las capas que son de medios idénticos consecutivos
%para disminuir el tiempo de calculo.
if Eliminacapas==1
    cont=1;
    daux(1)=d(1);
    deltaaux(1)=delta(1);
    Mediosaux(1)=Medios(1);
    for i=2:Nmed-2
        if(Medios(i+1)==Medios(i))
            daux(cont)=daux(cont)+d(i);
        else
            daux(cont+1)=d(i);
            deltaaux(cont+1)=delta(i);
            Mediosaux(cont+1)=Medios(i);
            cont=cont+1;
        end
    end
    %pone el ultimo medio
    Mediosaux(end+1:end+2)=Medios(end-1:end);
    deltaaux(end+1)=delta(end);
    
    %reemplazo por los vectores corregidos
    Medios=Mediosaux;
    d=daux;
    delta=deltaaux;
    Nmed=size(Medios,2);
end

%El programa da un mensaje de error y no realiza el cálculo si una rugosidad elegida es mayor que el
%espesor de una de las capas a un lado de la interfase en cuestión.
if rugosidad==1
    for l=1:length(d)
        if Ndelta*(delta(l)+delta(l+1))>= d(l);
            disp (['Error, la rugosidad elegida es mas grande que el elemento ', (num2str(l)), ' del vector de espesores.'])
            return
        end
    end
end

%Aqui podemos guardar o cargar el vector de espesores de la multicapa
save d.dat d -ascii
%d=load('d.dat');

%valores de teta
if Npt>1
    teta=tetai:(tetaf-tetai)/(Npt-1):tetaf;
else
    teta=tetai ;
end

%convierto a radianes
teta=teta/180*pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Realizo el cálculo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:Npl
    %creo el vector de indices de refracción para cada longitud de onda
    %elegida
    for k=1:Nmed
        Nl(k)=N(Medios(k),j);
    end
    if rugosidad_c==0
        if rugosidad==1
            [Nlrug, drug]=rugosear3(Nl,d,delta,Nepstein,Ndelta);
            %calculo el coeficiente de reflexion, transmitancia y la intensidad del
            %campo electrico al cuadrado respecto del incidente
            [R(j,:),T(j,:),I(:,:,j)]=campo5(Nlrug,drug,lambda(j),teta,p,ndiv);
            
        else
            %calculo el coeficiente de reflexion, transmitancia y la intensidad del
            %campo electrico al cuadrado respecto del incidente
            [R(j,:),T(j,:),I(:,:,j)]=campo5(Nl,d,lambda(j),teta,p,ndiv);
        end
    else
        R(j,:)=0;
        T(j,:)=0;
        for i=1:Nprom
            if rugosidad==1
                [Nlrug, drug]=rugosear3(Nl,d+randn(1,numel(d)).*delta_c,delta,Nepstein,Ndelta);
                %calculo el coeficiente de reflexion, transmitancia y la intensidad del
                %campo electrico al cuadrado respecto del incidente
                [Ri(j,:),Ti(j,:),I(:,:,j)]=campo5(Nlrug,drug,lambda(j),teta,p,ndiv);
                
            else
                %calculo el coeficiente de reflexion, transmitancia y la intensidad del
                %campo electrico al cuadrado respecto del incidente
                [Ri(j,:),Ti(j,:),I(:,:,j)]=campo5(Nl,d+randn(1,numel(d)).*delta_c,lambda(j),teta,p,ndiv);
            end
            R(j,:)=R(j,:)+Ri(j,:)/Nprom;
            T(j,:)=T(j,:)+Ti(j,:)/Nprom;
        end      
    end
    
end

%acomoda los resultados del campo electrico para que queden en forma de
%matriz.
I=squeeze(I);

%si se elijió graficar las multicapas se hace el grafico
if GraficaMulticapa==1
    %dibuja el perfil de indices de refraccion para la longitud de onda de
    %referencia
    for k=1:Nmed
        Nref(k)=N(Medios(k),Npl+1);
    end

    if rugosidad==1
        [Nref drug]=rugosear3(Nref,d,delta,Nepstein,Ndelta);
        Dibujaescalones2b(Nref,drug,Medios)
    else
        Dibujaescalones2b(Nref,d,Medios)
    end
end

%elimina el ultimo valor de lambda para que tenga la dimension de las
%reflectancias
lambda(Npl+1)=[];

%si tenemos un solo angulo grafico los epectros en funcion de lambda
if Npt==1
    figure(10)
    plot(lambda*1e9,R*100,'b.-')
    hold on
    plot(lambda*1e9,T*100,'.-g')
    plot(lambda*1e9,(1-R-T)*100,'.-r')
    xlabel('Longitud de onda [nm]')
    ylabel('Reflectancia [%]')
    legend('R ', 'T ','Abs')
    %legend('R')

salida=[lambda(:)*1e9 R T ];
save espectro.dat salida -ascii
toc
    
    if CalculodeCampo==1
        %dibujo el cociente cuadrado del campo electrico respecto del campo
        %incidente
        %calculo las posiciones donde se calculo el campo

        aux=d/ndiv;
        for k=1:ndiv-1
            aux=[aux;d/ndiv];
        end
        aux2=reshape(aux,size(aux,2)*ndiv,1);

        x=[0 cumsum(aux2')];
        
        figure
        contourf(lambda*1e9,x*1e6,log10(I))
        xlabel('Longitud de onda [nm]')
        ylabel('Posicion en la microcavidad [\mu m]')
        title('Logaritmo del cociente de intensidades')
        colorbar
        
        figure
        surf(lambda*1e9,x*1e6,I)
        xlabel('Longitud de onda [nm]')
        ylabel('Posicion en la microcavidad [\mu m]')
        zlabel('Cociente de intensidades')
        title('Campo eléctrico dentro de la multicapa')
        shading flat
        
        %dibuja el campo electrico en la microcavidad para la longitud de onda de
        %referencia superpuesta con la distribucion espacial del indice de
        %refraccion.
        %dibuja el perfil de indice de refraccion de las capas
        %si se pidio graficar las multicapas se hace el grafico
        
        %dibuja el perfil de indices de refraccion para la longitud de onda de
        %referencia
        for k=1:Nmed
            Nref(k)=N(Medios(k),Npl+1);
        end
        Dibujaescalones2b(Nref,d,Medios)
        
        hold on
        %calcula el campo para la longitud de onda de referencia
        [Rref,Tref,Iref(:,1)]=campo2(Nref,d,lambda0,teta,p,ndiv);
        plot(x*1e6,Iref)
    end
end

%si tenemos un solo lambda grafico en funcion de teta
if Npl==1
    teta=teta*180/pi;     %convierto a radianes
    figure
    plot(teta, R*100, '.', 'Color', [0.251, 0.878, 0.816])
    hold on
    grid on
    %plot(teta,T*100,'.g')
    %plot(teta,(1-R-T)*100,'.r')
    xlabel('Angulo de incidencia')
    ylabel('Reflectancia [%]')
    %legend('R ', 'T ','Abs')
    title("Curva teórica tipo TAMM")
    ylim([45 100])
    
    if CalculodeCampo==1
        %dibujo el cociente cuadrado del campo electrico respecto del campo
        %incidente
        %calculo las posiciones donde se calculo el campo
        
        aux=d/ndiv;
        for k=1:ndiv-1
            aux=[aux;d/ndiv];
        end
        aux2=reshape(aux,size(aux,2)*ndiv,1);
        
        x=[0 cumsum(aux2')];
        
        figure
        contourf(teta,x*1e6,log10(I'))
        xlabel('Angulo de incidencia')
        ylabel('Posicion en la microcavidad [\mu m]')
        title('Logaritmo del cociente de intensidades')
        colorbar
        
        figure
        surf(teta,x*1e6,I')
        xlabel('Angulo de incidencia')
        ylabel('Posicion en la microcavidad [\mu m]')
        shading flat
        
    end
end

%si tenemos mas de un teta y lambda dibujo en función de ambos
if and(Npl>1,Npt>1)
    teta=teta*180/pi;     %convierto a grados
    
    figure
    contourf(teta,lambda*1e9,R*100)
    ylabel('Longitud de onda [nm]')
    xlabel('Angulo de incidencia')
    title('Reflectancia')
    colorbar
    
    figure
    contourf(teta,lambda*1e9,T*100)
    ylabel('Longitud de onda [nm]')
    xlabel('Angulo de incidencia')
    title('Transmitancia')
    colorbar
    
    figure
    contourf(teta,lambda*1e9,(1-R-T)*100)
    ylabel('Longitud de onda [nm]')
    xlabel('Angulo de incidencia')
    title('Absorbancia')
    colorbar
end

end

%Definiciones de funciones usadas en el cálculo

function n = epstein(n1,n2,deltaep,z,z0)

% devuelve el perfil de refraccion de epstein para dos indices de
% refraccion limites n1 y n2 (n2>n1), el espesor de la homogeneidad d, la
% posicion z.
%
% Optics of thin film - An optical multilayer theory (Knittl, 1976).

n = sqrt(0.5*(n1.^2+n2.^2)+0.5*(n2.^2-n1.^2).*tanh((z-z0)./(deltaep/2)));

end

function n = epsteincomplejo(n1,n2,deltaep,z,z0)

nr= epstein(real(n1),real(n2),deltaep,z,z0);
ni= epstein(imag(n1),imag(n2),deltaep,z,z0);

n=nr+1i*ni;
end

function [n1,espesores]=rugosear3(n,esp,d,N,Nd)
%en esta version soporta valores de d=0
pos=[0 cumsum(esp)];

for j=1:length(d)
    if d(j)
        aux(j,:)=(-Nd*d(j):d(j)/N:Nd*d(j));
    end
end

n1=n(1);
x=[];

for i=1:length(d)
    if d(i)
        %generamos las posiciones
        aux1=pos(i)+aux(i,:);
        x=[x aux1];
        %generamos los valores de indice de refraccion
        n1=[n1 epsteincomplejo(n(i),n(i+1),d(i),aux(i,:),0)];
    else
         x=[x pos(i)];
         n1=[n1 n(i+1)];
    end
    
end

% figure
% plot(x,n1(2:end),'.-')
espesores=diff(x);
end

function [R,T,I]=campo5(N,d,lambda,teta,p,n)
%calculo de la reflectancia de una serie de capas dielectricas
%Tenemos Nc capas dielectricas de indice de refracción N(i). Este indice
%tiene partes real e imaginaria.
%en este caso para el calculo del campo se divide cada capa dielectrica en
%n partes iguales y se calcula el campo en cada una de ellas.
%teta puede ser vector. En ese
%caso R y T seran vectores e I sera una matriz. El valor de p indica si la
%polarización es s (p=0) o p (p=1). En esta version se aceptan valores
%intermedios de p.

if p>0
    %calculo para polarización p
    for m=1:length(teta)

        Nc=size(N,2);
        costeta=cos(teta(m));

        %calculo la admitancia del medio de llegada
        adm(1)=N(1)./cos(teta(m));

        %calculo el producto de todas las matrices internas correspondientes a las
        %Nc-2 capas dielectricas
        M=eye(2);
        Mc=eye(2);
        M1(1)=Mc(1,1);
        M2(1)=Mc(1,2);
        for j=2:Nc-1

            %calculo los angulos teta dentro de cada dielectrico usando la ley de
            %Snell
            %Para evitar problemas con el tratamiento de los angulos complejos
            %durante la reflexion total interna usamos el cálculo de los cosenos de
            %teta.

            costeta(j)=(1-(N(j-1)./N(j))^2*(1-costeta(j-1)^2))^0.5;

            %calculo los delta para cada capa (notar que son Nc-2)

            delta(j)=2*pi*N(j)*d(j-1)*costeta(j)/lambda;

            %calculo la admitancia de la capa correspondiente (la admitancia
            %depende de si la onda esta polarizada en forma p o en  forma s.

            adm(j)=N(j)./costeta(j);

            M=M*matriz(delta(j),adm(j));

            for k=1:n
                Mc=matrizC(delta(j)/n,adm(j))*Mc;
                M1(n*(j-2)+k+1)=Mc(1,1);
                M2(n*(j-2)+k+1)=Mc(1,2);
            end
        end

        %calculo la admitancia del ultimo medio
        %teta(Nc)=asin(sin(teta(Nc-1)*N(Nc-1)./N(Nc)));
        costeta(Nc)=(1-(N(Nc-1)./N(Nc))^2*(1-costeta(Nc-1)^2))^0.5;

        adm(Nc)=N(Nc)./costeta(Nc);

        %construida la matriz calculo la admitancia equivalente de la
        %multicapa
        aux=M*[1 adm(Nc)]';

        Y=aux(2)./aux(1);

        %calculo el coeficiente de reflexion
        r=(adm(1)-Y)./(adm(1)+Y);

        %calculo el coeficiente de transmicion
        C=aux(2);
        B=aux(1);
        Tp(m)=real(4*adm(1)*adm(Nc)/(abs(adm(1)*B+C)).^2);

        %calculo la reflectancia
        Rp(m)=abs(r).^2;

        %calculo la relacion de intensidades del campo electrico al cuadrado
        %(normalizado con el campo incidente)
        Ip(m,:)=abs(M1*(1+r)+M2*adm(1)*(1-r)).^2;
    end
else
    Tp=0;
    Rp=0;
    Ip=0;
end


if p<1
    %calculo para polarización s
    for m=1:length(teta)

        Nc=size(N,2);
        costeta=cos(teta(m));

        %calculo la admitancia del medio de llegada
        adm(1)=N(1)*cos(teta(m));

        %calculo el producto de todas las matrices internas correspondientes a las
        %Nc-2 capas dielectricas
        M=eye(2);
        Mc=eye(2);
        M1(1)=Mc(1,1);
        M2(1)=Mc(1,2);
        for j=2:Nc-1

            %calculo los angulos teta dentro de cada dielectrico usando la ley de
            %Snell
            %Para evitar problemas con el tratamiento de los angulos complejos
            %durante la reflexion total interna usamos el cálculo de los cosenos de
            %teta.

            costeta(j)=(1-(N(j-1)./N(j))^2*(1-costeta(j-1)^2))^0.5;

            %calculo los delta para cada capa (notar que son Nc-2)

            delta(j)=2*pi*N(j)*d(j-1)*costeta(j)/lambda;

            %calculo la admitancia de la capa correspondiente (la admitancia
            %depende de si la onda esta polarizada en forma p o en  forma s.

            adm(j)=N(j)*costeta(j);

            M=M*matriz(delta(j),adm(j));

            for k=1:n
                Mc=matrizC(delta(j)/n,adm(j))*Mc;
                M1(n*(j-2)+k+1)=Mc(1,1);
                M2(n*(j-2)+k+1)=Mc(1,2);
            end
        end

        %calculo la admitancia del ultimo medio
        %teta(Nc)=asin(sin(teta(Nc-1)*N(Nc-1)./N(Nc)));
        costeta(Nc)=(1-(N(Nc-1)./N(Nc))^2*(1-costeta(Nc-1)^2))^0.5;

        adm(Nc)=N(Nc)*costeta(Nc);

        %construida la matriz calculo la admitancia equivalente de la
        %multicapa
        aux=M*[1 adm(Nc)]';

        Y=aux(2)./aux(1);

        %calculo el coeficiente de reflexion
        r=(adm(1)-Y)./(adm(1)+Y);

        %calculo el coeficiente de transmicion
        C=aux(2);
        B=aux(1);
        Ts(m)=real(4*adm(1)*adm(Nc)/(abs(adm(1)*B+C)).^2);

        %calculo la reflectancia
        Rs(m)=abs(r).^2;

        %calculo la relacion de intensidades del campo electrico al cuadrado
        %(normalizado con el campo incidente)
        Is(m,:)=abs(M1*(1+r)+M2*adm(1)*(1-r)).^2;
    end
else
    Ts=0;
    Rs=0;
    Is=0;
end

R=p*Rp+(1-p)*Rs;
T=p*Tp+(1-p)*Ts;
I=p*Ip+(1-p)*Is;
%deberia chequearse que la relacion de los campos electricos al cuadrado
%se suma de esta manera.
end

function [R,T,I]=campo2(N,d,lambda,teta,p,n)
%calculo de la reflectancia de una serie de capas dielectricas
%en esta version se solucionó el problema de la reflexion total de la
%version 1. Tenemos Nc capas dielectricas de indice de refracción N(i). Este indice
%tiene partes real e imaginaria.
%en este caso para el calculo del campo se divide cada capa dielectrica en
%n partes iguales y se calcula el campo en cada una de ellas.

Nc=size(N,2);
costeta=cos(teta);

%calculo la admitancia del medio de llegada
if p==1
    adm(1)=N(1)./cos(teta(1));
else
    adm(1)=N(1)*cos(teta(1));
end


%calculo el producto de todas las matrices internas correspondientes a las
%Nc-2 capas dielectricas
M=eye(2);
Mc=eye(2);
M1(1)=Mc(1,1);
M2(1)=Mc(1,2);
for j=2:Nc-1

    %calculo los angulos teta dentro de cada dielectrico usando la ley de
    %Snell
    %Para evitar problemas con el tratamiento de los angulos complejos
    %durante la reflexion total interna usamos el cálculo de los cosenos de
    %teta.

    costeta(j)=(1-(N(j-1)./N(j))^2*(1-costeta(j-1)^2))^0.5;

    %calculo los delta para cada capa (notar que son Nc-2)

    delta(j)=2*pi*N(j)*d(j-1)*costeta(j)/lambda;

    %calculo la admitancia de la capa correspondiente (la admitancia
    %depende de si la onda esta polarizada en forma p o en  forma s.

    if p==1
        adm(j)=N(j)./costeta(j);
    else
        adm(j)=N(j)*costeta(j);
    end

    M=M*matriz(delta(j),adm(j));

    for k=1:n
        Mc=matrizC(delta(j)/n,adm(j))*Mc;
        M1(n*(j-2)+k+1)=Mc(1,1);
        M2(n*(j-2)+k+1)=Mc(1,2);
    end
end

%calculo la admitancia del ultimo medio
%teta(Nc)=asin(sin(teta(Nc-1)*N(Nc-1)./N(Nc)));
costeta(Nc)=(1-(N(Nc-1)./N(Nc))^2*(1-costeta(Nc-1)^2))^0.5;

if p==1
    adm(Nc)=N(Nc)./costeta(Nc);
else
    adm(Nc)=N(Nc)*costeta(Nc);
end


%construida la matriz calculo calculo la admitancia equivalente de la
%multicapa
aux=M*[1 adm(Nc)]';

Y=aux(2)./aux(1);

%calculo el coeficiente de reflexion
r=(adm(1)-Y)./(adm(1)+Y);

%calculo el coeficiente de transmicion
C=aux(2);
B=aux(1);
T=real(4*adm(1)*adm(Nc)/(abs(adm(1)*B+C)).^2);

%calculo la reflectancia
R=abs(r).^2;

%calculo la relacion de intensidades del campo electrico al cuadrado
%(normalizado con el campo incidente)
I=abs(M1*(1+r)+M2*adm(1)*(1-r)).^2;
end

function M=matriz(delta,adm)
%genero la matriz de  la capa correspondiente

M(1,1)=cos(delta);
M(1,2)=-1i/adm*sin(delta);
M(2,1)=-1i*adm*sin(delta);
M(2,2)=cos(delta);
end

function M=matrizC(delta,adm)
%genero la matriz de  la capa correspondiente

M(1,1)=cos(delta);
M(1,2)=1i/adm*sin(delta);
M(2,1)=1i*adm*sin(delta);
M(2,2)=cos(delta);
end

function Dibujaescalones2b(N,d,Medios)
%dibuja el perfil de indice de refraccion elegido apartir de los valores de
%desplazamiento e indice

aux=[-0.1*sum(d) cumsum([0 d ])]; %le asigno un 0.1 del espesor total de la multicapa al medio de la izquierda
x=reshape(cat(1,aux,aux),size(aux,2)*2,1);
y=reshape(cat(1,N,N),size(N,2)*2,1);

x(1)=[];
x(size(x,1)+1)= x(size(x,1))+ 0.1*sum(d); %le asigno un 0.1 del espesor total de la multicapa al medio de la derecha
figure
% plot(x*1e6,real(y))
% xlabel('distancia [\mum]')
% ylabel('Indice de refraccion real')
%axis tight
hold on

au=zeros(size(y,1),size(y,1)/2);
for j=1:size(y,1)/2
    au(j*2-1,j)=real(y(j*2-1));
    au(j*2,j)=real(y(j*2));
end

%creo el mapa de colores segun el medio utilizado (se pueden utilizar hasta
%64 medios diferentes)
% aux2=colormap;
% for j=1:size(Medios,2)
%     map(j,:)=aux2(Medios(j)*floor(64/max(Medios)),:);
% end
% colormap(map);

%leyenda=unique(map,'rows');
area(x*1e6,au,'EdgeColor','none');
%aux3=[1 2 3 4 21];
%h=area(x*1e6,au(:,aux3));

%legend(h,'Vidrio','SnO2','Si 82%','Si 55%','Aluminio')

end

%%%%%%%%%% Definiciones de funciones de indice de refraccion %%%%%%%%%%%%%%%%%%%%

function N=dummy(lambda,n)
%función que devuelve un indice de refracción constante
N=lambda*0+n;
end

function N=indiceaguaLyT(lambda,temperatura)
%funcion que devuelve el indice de refraccion  del  agua
%esta funcion toma la temperatura como argumento.
%La temperatura esta en grados centrigrados y la longitud de onda en metros

lambdas=[226.50  361.05  404.41  589.00  632.80 1013.98 ]*1e-9;
temperaturas=[0    10    20    30    40    50    60    70    80    90   100];

aux=    [1.3945    1.3490    1.3442    1.3343    1.3331    1.3261
    1.3942    1.3487    1.3439    1.3341    1.3328    1.3259
    1.3934    1.3480    1.3432    1.3334    1.3321    1.3252
    1.3921    1.3468    1.3421    1.3323    1.3311    1.3242
    1.3905    1.3454    1.3406    1.3310    1.3297    1.3230
    1.3885    1.3437    1.3390    1.3294    1.3281    1.3215
    1.3864    1.3418    1.3371    1.3276    1.3264    1.3197
    1.3840    1.3397    1.3351    1.3256    1.3244    1.3178
    1.3813    1.3375    1.3328    1.3234    1.3222    1.3158
    1.3785    1.3350    1.3304    1.3211    1.3199    1.3135
    1.3755    1.3324    1.3278    1.3186    1.3174    1.3111];

N=interp2(lambdas,temperaturas,aux,lambda,temperatura,'cubic');
end

function N=indiceagua(lambda)
%funcion que devuelve el indice de refraccion  del  agua
%Tomado de  Segelstein, D., 1981: "The Complex Refractive Index of Water",
% M.S. Thesis, University of Missouri--Kansas City
%Wavelength(um)	Real Index	Imag Indexs
datos=load('agua.txt');
Nr=interp1(datos(:,1)*1e-6,datos(:,2),lambda);
Ni=interp1(datos(:,1)*1e-6,datos(:,3),lambda);
N=Nr+1i*Ni;
end

function N=indicealcohol(lambda)
%funcion que devuelve el indice de refraccion  del alcohol ethilico
%ver http://refractiveindex.info/

lambda=lambda*1e6; %las formulas estan en micrómetros

C1 = 1.35265;
C2 = 0.00306;
C3 = 0.00002;

N= C1 + C2./(lambda.^2) + C3./(lambda.^4);

end

function N=indiceAIP(lambda)
%funcion que devuelve el indice de refraccion  del alcohol isopropilico
%Obtenido de https://refractiveindex.info/?shelf=organic&book=propanol&page=Kozma
%I. Z. Kozma, P. Krok, and E. Riedle. Direct measurement of the group-velocity mismatch and derivation of the refractive-index dispersion for a variety of solvents in the ultraviolet, J. Opt. Soc. Am. B 22, 1479-1485 (2005)
lambda=lambda*1e6;

N=1.36485+4.29404081e-3.*lambda.^-2-6.4823380e-5.*lambda.^-4+3.41833e-6.*lambda.^-6;
end

function n=ZnO(lambda)
%indice de refraccion del oxido de Zn obtenido de
%http://refractiveindex.info/?group=CRYSTALS&material=ZnO
%Es el rayo ordinario y la longitud de onda en metros

C1 = 2.81418;
C2 = 0.87968;
C3 = 0.3042;
C4 = 0.00711;
lambda=lambda*1e6;

n = sqrt(C1 + C2*lambda.^2./(lambda.^2-C3^2) - C4*lambda.^2);

%el rayo extraordinario tiene otros parametros
%C1 = 2.80333; C2 = 0.94470; C3 = 0.3004; C4 = 0.00714
end

function N=Al2O3(lambda)
%indice de refraccion del aluminio obtenido de
%http://refractiveindex.info/?group=CRYSTALS&material=Al2O3

longdaexp=[2.48	2.066	1.771	1.55	1.378	1.24	1.127	1.033	0.9537	0.8856	0.8266	0.7749	0.7293	0.6888	0.6525	0.6199	0.5904	0.5636	0.5391	0.5166	0.4959	0.4769	0.4592	0.4428	0.4275	0.4133	0.3999	0.3875	0.3757	0.3647	0.3542	0.3444	0.3351	0.3263	0.3179	0.31	0.3024	0.2952	0.2883	0.2818	0.2755	0.2695	0.2638	0.2583	0.253	0.248	0.2431	0.2384	0.2339	0.2296	0.2254	0.2214	0.2175	0.2138	0.2101	0.2066];
n=[1.726	1.736	1.742	1.746	1.749	1.751	1.7526	1.7542	1.7558	1.7574	1.759	1.7606	1.7622	1.7638	1.7654	1.767	1.76857	1.77014	1.77174	1.77335	1.775	1.77666	1.77837	1.78015	1.78202	1.784	1.78619	1.7885	1.7909	1.79341	1.796	1.79861	1.8013	1.8041	1.80699	1.81	1.81358	1.81717	1.82067	1.82398	1.827	1.8291	1.83083	1.83221	1.83326	1.834	1.83434	1.83445	1.83439	1.83422	1.834	1.83378	1.83361	1.83355	1.83366	1.834];
aux=interp1(longdaexp*1e-6,n,lambda);

N=aux+1i*0;


end

function N=Al2O3amorfilm(lambda)
%indice de refraccion de la alumina obtenida por deposicion 

% REFERENCES: "R. Boidin, T. Halenkovi?, V. Nazabal, L. Bene, P. N?mec. Pulsed laser deposited alumina thin films, <a href=\"http://dx.doi.org/10.1016/j.ceramint.2015.09.048\"><i>Ceramics International</i> <b>42</b>, 1177-1182 (2016)</a> (Numerical data kindly provided by Tomá Halenkovi?)"
% COMMENTS: "Thin films of amorphous alumina were fabricated using pulsed laser deposition in vacuum.<br>Influence of deposition parameters, specifically deposition time under vacuum and partial pressure of argon, is discussed in the original publication."
% DATA:

longdaexp=[0.3	0.32	0.34	0.36	0.38	0.4	0.42	0.44	0.46	0.48	0.5	0.52	0.54	0.56	0.58	0.6	0.62	0.64	0.66	0.68	0.7	0.72	0.74	0.76	0.78	0.8	0.82	0.84	0.86	0.88	0.9	0.92	0.94	0.96	0.98	1	1.02	1.04	1.06	1.08	1.1	1.12	1.14	1.16	1.18	1.2	1.22	1.24	1.26	1.28	1.3	1.32	1.34	1.36	1.38	1.4	1.42	1.44	1.46	1.48	1.5	1.52	1.54	1.56	1.58	1.6	1.62	1.64	1.66	1.68	1.7	1.72	1.74	1.76	1.78	1.8	1.82	1.84	1.86	1.88	1.9	1.92	1.94	1.96	1.98	2	2.02	2.04	2.06	2.08	2.1	2.12	2.14	2.16	2.18	2.2	2.22	2.24	2.26	2.28	2.3	2.3188	2.3398	2.3568	2.3784	2.396	2.4183	2.4365	2.4596	2.4784	2.4975	2.5169	2.5366	2.5567	2.577	2.5976	2.6186	2.64	2.6562	2.6782	2.6949	2.7175	2.7347	2.7579	2.7756	2.7996	2.8179	2.8364	2.8551	2.8741	2.8998	2.9194	2.9393	2.9594	2.9798	2.9936	3.0145	3.0357	3.0571	3.0789	3.0936	3.1159	3.1386	3.1538	3.177	3.1927	3.2164	3.2325	3.2569	3.2733	3.2983	3.3152	3.3322	3.3581	3.3756	3.3933	3.4111	3.4383	3.4566	3.4751	3.4939	3.5128	3.532	3.5513	3.5709	3.5907	3.6107	3.6309	3.6513	3.672	3.693	3.7141	3.7355	3.7572	3.7791	3.7901	3.8124	3.835	3.8578	3.8693	3.8926	3.9161	3.9399	3.9519	3.9762	3.9884	4.0131	4.0381	4.0507	4.0762	4.089	4.115	4.1281	4.1546	4.1679	4.1949	4.2085	4.236	4.2499	4.278	4.2921	4.3064	4.3352	4.3498	4.3791	4.394	4.4089	4.4391	4.4544	4.4697	4.4852	4.5165	4.5323	4.5482	4.5642	4.5965	4.6129	4.6294	4.646	4.6795	4.6965	4.7136	4.7308	4.7481	4.7655	4.7831	4.8187	4.8367	4.8548	4.873	4.8914	4.91	4.9286	4.9474	4.9664	4.9855	5.0047	5.0241	5.0437	5.0634	5.0832	5.1033	5.1234	5.1438	5.1643	5.1849	5.2057	5.2267	5.2479	5.2692	5.2907	5.3124	5.3343	5.3563	5.3785	5.4009	5.4235	5.4463	5.4693	5.4925	5.5159	5.5394	5.5632	5.5872	5.6114	5.6358	5.6604	5.6852	5.7103	5.7355	5.761	5.7867	5.8127	5.8389	5.8653	5.8919	5.9188	5.946	5.9734	6.0011	6.029	6.0571	6.0856	6.1143	6.1433	6.1725	6.202	6.2319	6.262	6.2924	6.3231	6.3541	6.3854	6.417	6.4489	6.4811	6.5137	6.5466	6.5798	6.6134	6.6473	6.6816	6.7162	6.7512	6.7865	6.8223	6.8583	6.8948	6.9317	6.969	7.0066	7.0447	7.0832	7.1221	7.1615	7.2013	7.2415	7.2822	7.3233	7.3649	7.407	7.4496	7.4926	7.5362	7.5803	7.6249	7.67	7.7156	7.7618	7.8086	7.8559	7.9038	7.9523	8.0014	8.0511	8.1014	8.1524	8.204	8.2562	8.3091	8.3628	8.4171	8.4721	8.5278	8.5843	8.6415	8.6995	8.7583	8.8179	8.8783	8.9395	9.0016	9.0645	9.1284	9.1931	9.2588	9.3254	9.393	9.4615	9.5311	9.6017	9.6733	9.7461	9.8199	9.8949	9.971	10.048	10.127	10.207	10.288	10.37	10.453	10.538	10.625	10.713	10.802	10.893	10.985	11.079	11.174	11.272	11.37	11.471	11.573	11.678	11.784	11.892	12.002	12.114	12.229	12.345	12.464	12.585	12.708	12.834	12.962	13.093	13.227	13.363	13.502	13.645	13.79	13.938	14.089	14.244	14.403	14.564	14.73	14.899	15.072	15.25	15.431	15.617	15.808	16.003	16.203	16.408	16.618	16.834	17.056	17.283	17.517	17.757	18.003];
n=[1.73756	1.72717	1.71886	1.7121	1.70652	1.70185	1.69791	1.69454	1.69164	1.68912	1.68691	1.68496	1.68324	1.68169	1.68031	1.67906	1.67792	1.67689	1.67594	1.67507	1.67427	1.67353	1.67283	1.67219	1.67158	1.67101	1.67047	1.66996	1.66948	1.66902	1.66858	1.66816	1.66776	1.66737	1.66699	1.66663	1.66628	1.66594	1.66561	1.66528	1.66496	1.66465	1.66435	1.66405	1.66376	1.66347	1.66318	1.6629	1.66262	1.66234	1.66207	1.6618	1.66152	1.66125	1.66099	1.66072	1.66045	1.66019	1.65992	1.65965	1.65939	1.65912	1.65886	1.65859	1.65832	1.65805	1.65778	1.65751	1.65724	1.65697	1.6567	1.65642	1.65615	1.65587	1.65559	1.65531	1.65503	1.65475	1.65446	1.65418	1.65389	1.6536	1.6533	1.65301	1.65271	1.65241	1.65211	1.65181	1.6515	1.6512	1.65089	1.65057	1.65026	1.64994	1.64962	1.6493	1.64898	1.64865	1.64832	1.64799	1.64765	1.64733	1.64698	1.64669	1.64632	1.64601	1.64562	1.6453	1.64489	1.64455	1.6442	1.64385	1.64349	1.64312	1.64274	1.64235	1.64195	1.64155	1.64123	1.64081	1.64048	1.64004	1.6397	1.63923	1.63887	1.63839	1.63801	1.63763	1.63724	1.63685	1.6363	1.63589	1.63546	1.63503	1.63458	1.63428	1.63382	1.63335	1.63287	1.63238	1.63205	1.63153	1.63101	1.63066	1.63012	1.62975	1.62918	1.6288	1.62821	1.62781	1.6272	1.62679	1.62636	1.62572	1.62528	1.62483	1.62437	1.62368	1.6232	1.62272	1.62223	1.62172	1.62122	1.6207	1.62017	1.61963	1.61908	1.61853	1.61796	1.61738	1.61679	1.61619	1.61558	1.61496	1.61432	1.614	1.61335	1.61268	1.612	1.61165	1.61095	1.61024	1.60951	1.60914	1.60839	1.608	1.60723	1.60644	1.60603	1.60522	1.6048	1.60396	1.60353	1.60266	1.60221	1.60131	1.60086	1.59992	1.59945	1.59849	1.598	1.5975	1.59649	1.59598	1.59494	1.59441	1.59387	1.59278	1.59222	1.59166	1.59109	1.58993	1.58933	1.58874	1.58813	1.5869	1.58627	1.58563	1.58498	1.58367	1.583	1.58232	1.58163	1.58093	1.58023	1.57951	1.57805	1.57731	1.57655	1.57579	1.57501	1.57423	1.57343	1.57262	1.5718	1.57097	1.57013	1.56928	1.56841	1.56753	1.56664	1.56573	1.56482	1.56388	1.56294	1.56198	1.561	1.56002	1.55901	1.55799	1.55696	1.55591	1.55484	1.55375	1.55265	1.55153	1.5504	1.54924	1.54807	1.54687	1.54566	1.54443	1.54317	1.5419	1.5406	1.53928	1.53794	1.53657	1.53518	1.53377	1.53233	1.53087	1.52937	1.52786	1.52631	1.52473	1.52313	1.5215	1.51983	1.51813	1.5164	1.51464	1.51284	1.511	1.50913	1.50722	1.50527	1.50329	1.50126	1.49918	1.49707	1.4949	1.4927	1.49044	1.48813	1.48577	1.48336	1.48089	1.47837	1.47579	1.47315	1.47044	1.46767	1.46483	1.46192	1.45894	1.45589	1.45276	1.44954	1.44625	1.44286	1.43939	1.43582	1.43216	1.42839	1.42452	1.42054	1.41644	1.41223	1.40789	1.40342	1.39882	1.39408	1.38919	1.38414	1.37893	1.37356	1.36801	1.36227	1.35634	1.35021	1.34386	1.33728	1.33047	1.3234	1.31607	1.30846	1.30055	1.29233	1.28377	1.27485	1.26556	1.25586	1.24573	1.23514	1.22405	1.21243	1.20023	1.18742	1.17394	1.15975	1.14477	1.12895	1.11222	1.0945	1.0757	1.05574	1.03453	1.01195	0.98791	0.9623	0.93501	0.90598	0.87515	0.84254	0.80828	0.77268	0.73636	0.70035	0.66616	0.63567	0.6108	0.59289	0.58242	0.57913	0.58231	0.59116	0.60493	0.62301	0.64491	0.67024	0.6987	0.73004	0.76407	0.8006	0.83946	0.88053	0.92365	0.96867	1.01545	1.06385	1.11369	1.16482	1.21708	1.27028	1.32423	1.37875	1.43366	1.48873	1.54378	1.5986	1.65299	1.70675	1.75969	1.81159	1.86229	1.91161	1.95938	2.00542	2.04962	2.09184	2.13195	2.16987	2.20551	2.23881	2.26973	2.29825	2.32436	2.34808	2.36943	2.38848	2.4053	2.41997	2.43259];
aux=interp1(longdaexp*1e-6,n,lambda,'linear','extrap');

N=aux+1i*0;

end

function N=Ag(lambda)
%Indice de refraccion complejo de la plata metalica

%REFERENCES: "P. B. Johnson and R. W. Christy. Optical Constants of the Noble Metals, <a href=\"https://doi.org/10.1103/PhysRevB.6.4370\"><i>Phys. Rev. B</i> <b>6</b>, 4370-4379 (1972)</a>"
%COMMENTS: "Room temperature"

longdaexp=[0.1879	0.1916	0.1953	0.1993	0.2033	0.2073	0.2119	0.2164	0.2214	0.2262	0.2313	0.2371	0.2426	0.249	0.2551	0.2616	0.2689	0.2761	0.2844	0.2924	0.3009	0.3107	0.3204	0.3315	0.3425	0.3542	0.3679	0.3815	0.3974	0.4133	0.4305	0.4509	0.4714	0.4959	0.5209	0.5486	0.5821	0.6168	0.6595	0.7045	0.756	0.8211	0.892	0.984	1.088	1.216	1.393	1.61	1.937];
n=[1.07	1.1	1.12	1.14	1.15	1.18	1.2	1.22	1.25	1.26	1.28	1.28	1.3	1.31	1.33	1.35	1.38	1.41	1.41	1.39	1.34	1.13	0.81	0.17	0.14	0.1	0.07	0.05	0.05	0.05	0.04	0.04	0.05	0.05	0.05	0.06	0.05	0.06	0.05	0.04	0.03	0.04	0.04	0.04	0.04	0.09	0.13	0.15	0.24];
k=[1.212	1.232	1.255	1.277	1.296	1.312	1.325	1.336	1.342	1.344	1.357	1.367	1.378	1.389	1.393	1.387	1.372	1.331	1.264	1.161	0.964	0.616	0.392	0.829	1.142	1.419	1.657	1.864	2.07	2.275	2.462	2.657	2.869	3.093	3.324	3.586	3.858	4.152	4.483	4.838	5.242	5.727	6.312	6.992	7.795	8.828	10.1	11.85	14.08];
auxn=interp1(longdaexp*1e-6,n,lambda,'linear','extrap');
auxk=interp1(longdaexp*1e-6,k,lambda,'linear','extrap');

N=auxn+1i*auxk;

end

function N=Au(lambda)
%Indice de refraccion complejo del Oro metalico

%REFERENCES: "P. B. Johnson and R. W. Christy. Optical Constants of the Noble Metals, <a href=\"https://doi.org/10.1103/PhysRevB.6.4370\"><i>Phys. Rev. B</i> <b>6</b>, 4370-4379 (1972)</a>"
%COMMENTS: "Room temperature"

longdaexp=[0.1879	0.1916	0.1953	0.1993	0.2033	0.2073	0.2119	0.2164	0.2214	0.2262	0.2313	0.2371	0.2426	0.249	0.2551	0.2616	0.2689	0.2761	0.2844	0.2924	0.3009	0.3107	0.3204	0.3315	0.3425	0.3542	0.3679	0.3815	0.3974	0.4133	0.4305	0.4509	0.4714	0.4959	0.5209	0.5486	0.5821	0.6168	0.6595	0.7045	0.756	0.8211	0.892	0.984	1.088	1.216	1.393	1.61	1.937];
n=[1.28	1.32	1.34	1.33	1.33	1.3	1.3	1.3	1.3	1.31	1.3	1.32	1.32	1.33	1.33	1.35	1.38	1.43	1.47	1.49	1.53	1.53	1.54	1.48	1.48	1.5	1.48	1.46	1.47	1.46	1.45	1.38	1.31	1.04	0.62	0.43	0.29	0.21	0.14	0.13	0.14	0.16	0.17	0.22	0.27	0.35	0.43	0.56	0.92];
k=[1.188	1.203	1.226	1.251	1.277	1.304	1.35	1.387	1.427	1.46	1.497	1.536	1.577	1.631	1.688	1.749	1.803	1.847	1.869	1.878	1.889	1.893	1.898	1.883	1.871	1.866	1.895	1.933	1.952	1.958	1.948	1.914	1.849	1.833	2.081	2.455	2.863	3.272	3.697	4.103	4.542	5.083	5.663	6.35	7.15	8.145	9.519	11.21	13.78];
auxn=interp1(longdaexp*1e-6,n,lambda,'linear','extrap');
auxk=interp1(longdaexp*1e-6,k,lambda,'linear','extrap');

N=auxn+1i*auxk;

end

function N=Pt(lambda)
%Indice de refraccion complejo del Platino metalico

%REFERENCES: "W. S. M. Werner, K. Glantschnig, C. Ambrosch-Draxl. Optical constants and inelastic electron-scattering data for 17 elemental metals, <a href=\"https://doi.org/10.1063/1.3243762\"><i>J. Phys Chem Ref. Data</i> <b>38</b>, 1013-1092 (2009)</a>"
%COMMENTS: "Experimental data: Derived from reflection electron energy-loss spectroscopy (REELS) spectra."

longdaexp=[0.01759	0.01771	0.01784	0.01797	0.0181	0.01823	0.01837	0.01851	0.01864	0.01879	0.01893	0.01907	0.01922	0.01937	0.01953	0.01968	0.01984	0.02	0.02016	0.02032	0.02049	0.02066	0.02084	0.02101	0.02119	0.02138	0.02156	0.02175	0.02194	0.02214	0.02234	0.02254	0.02275	0.02296	0.02318	0.02339	0.02362	0.02384	0.02407	0.02431	0.02455	0.0248	0.02505	0.0253	0.02556	0.02583	0.0261	0.02638	0.02666	0.02695	0.02725	0.02755	0.02786	0.02818	0.0285	0.02883	0.02917	0.02952	0.02988	0.03024	0.03061	0.031	0.03139	0.03179	0.0322	0.03263	0.03306	0.03351	0.03397	0.03444	0.03492	0.03542	0.03594	0.03647	0.03701	0.03757	0.03815	0.03875	0.03936	0.04	0.04065	0.04133	0.04203	0.04275	0.0435	0.04428	0.04509	0.04592	0.04679	0.04769	0.04862	0.04959	0.05061	0.05166	0.05276	0.05391	0.0551	0.05636	0.05767	0.05904	0.06048	0.06199	0.06358	0.06525	0.06702	0.06888	0.07085	0.07293	0.07514	0.07749	0.07999	0.08266	0.08551	0.08856	0.09184	0.09537	0.09919	0.10332	0.10781	0.11271	0.11808	0.12398	0.13051	0.13776	0.14586	0.15498	0.16531	0.17712	0.19075	0.20664	0.22543	0.24797	0.26102	0.27552	0.29173	0.30996	0.33063	0.35424	0.38149	0.41328	0.45085	0.49594	0.55104	0.61992	0.70848	0.82656	0.99187	1.23984	1.65312	2.47968];
n=[0.8718	0.8709	0.8699	0.8695	0.8685	0.8675	0.8661	0.8651	0.8636	0.8621	0.8606	0.8591	0.857	0.855	0.8523	0.8503	0.8473	0.8441	0.8411	0.8376	0.8341	0.8302	0.8258	0.8215	0.8171	0.8128	0.8091	0.8061	0.8049	0.8071	0.813	0.8241	0.8408	0.8623	0.8844	0.9034	0.915	0.9182	0.9124	0.8974	0.8681	0.8233	1.0933	1.0283	0.9782	0.9517	0.9333	0.9187	0.9068	0.896	0.8858	0.8762	0.8677	0.8588	0.8504	0.8427	0.8346	0.8272	0.82	0.8136	0.8086	0.8044	0.8006	0.7979	0.7954	0.7913	0.7858	0.7787	0.7698	0.7603	0.7504	0.7393	0.7283	0.7174	0.7068	0.6958	0.6853	0.6766	0.6689	0.6642	0.664	0.6706	0.6866	0.7135	0.7468	0.777	0.7921	0.7903	0.7764	0.7571	0.738	0.7206	0.7069	0.6965	0.6903	0.6889	0.6925	0.7025	0.7197	0.7457	0.7829	0.8332	0.8981	0.9756	1.0582	1.132	1.1831	1.2052	1.2027	1.1851	1.1619	1.1398	1.1224	1.1124	1.11	1.1157	1.1288	1.149	1.1754	1.2072	1.2427	1.2857	1.371	1.4624	1.4746	1.4658	1.459	1.433	1.3774	1.2876	1.1748	1.1058	1.1449	1.2192	1.1883	1.1505	1.304	1.5367	1.2882	0.8675	0.6273	0.5124	0.4643	0.4611	0.5013	0.5979	0.787	1.1594	1.9756	4.2255];
k=[0.1182	0.12	0.1213	0.1225	0.1238	0.1251	0.127	0.1283	0.1297	0.1311	0.1325	0.1344	0.1359	0.138	0.1396	0.1417	0.1446	0.1469	0.1498	0.1534	0.1571	0.162	0.1671	0.1729	0.1805	0.1889	0.199	0.2115	0.2255	0.2416	0.2589	0.276	0.2896	0.2975	0.2968	0.2884	0.2743	0.2592	0.246	0.2396	0.2482	0.3298	0.3787	0.2057	0.184	0.1781	0.1763	0.1763	0.177	0.1786	0.1806	0.1837	0.1867	0.191	0.1952	0.2005	0.2061	0.2128	0.2201	0.228	0.2362	0.2449	0.2529	0.2601	0.2659	0.2704	0.2749	0.28	0.2858	0.2933	0.3018	0.3124	0.3247	0.3387	0.3544	0.373	0.3933	0.4168	0.4432	0.4735	0.5068	0.5428	0.5782	0.6075	0.6227	0.621	0.6078	0.5947	0.5906	0.5977	0.6145	0.639	0.6698	0.7057	0.746	0.7903	0.8382	0.8897	0.9434	0.999	1.055	1.1078	1.1518	1.1793	1.1827	1.1586	1.1153	1.0675	1.0293	1.0087	1.0065	1.0208	1.0478	1.0842	1.1275	1.1751	1.2252	1.2764	1.328	1.3792	1.4315	1.4953	1.5237	1.5785	1.5004	1.522	1.5459	1.5768	1.6288	1.7257	1.9066	2.2174	2.3998	2.5123	2.6012	2.8379	3.1046	3.0624	2.908	3.2129	3.7605	4.3965	5.121	5.9757	7.0273	8.3824	10.2266	12.9215	17.2802	25.5617];

auxn=interp1(longdaexp*1e-6,n,lambda,'linear','extrap');
auxk=interp1(longdaexp*1e-6,k,lambda,'linear','extrap');

N=auxn+1i*auxk;

end

function N=PZT(lambda)
%Indice de refraccion complejo del PZT obtenido por sol gel

%REFERENCES: Digitalizado de : J Mater Sci: Mater Electron
%DOI 10.1007/s10854-015-2859-9

longdaexp_n=[240.51034	245.64861	252.05707	263.56549	274.44765	283.42378	289.19346	294.32804	302.00744	306.49219	311.59583	316.05847	319.87938	324.33171	328.14673	331.95291	341.4713	354.14972	369.36264	393.48169	416.35191	437.32935	465.31058	501.56799	541.01343	594.46789	635.83779	679.12108	738.95307	784.78118	822.96895	840.79079];
n=[23.82387	24.46213	25.05918	25.7591	26.60318	27.49877	28.06495	28.65173	29.21787	29.61934	29.77368	29.86626	29.88678	29.83522	29.77337	29.58797	29.16565	28.4241	27.51775	26.52888	25.87979	25.46756	25.07578	24.69415	24.38453	24.13642	24.03268	23.98039	23.8969	23.82397	23.73059	23.69937];

longdaexp_k=[242.86478	252.42296 262.61032		282.34998	291.25165	296.3302	302.03424	307.10394	315.34393 319.1353	326.1051	331.17554	338.15567	350.86127	363.58017	376.31235	389.68628	401.7915	417.72278	440.02646	460.41965	480.8114	505.02698	526.05896	573.2207	605.7236	633.1279	690.48224	720.4336	746.56116	773.9648	794.99603	821.7631	849.1682];
k=[8.69091	8.65983	8.45929	8.10443	7.67286	7.25679	6.59424	5.99332 5.05342	4.37552	3.52808	2.94256	2.31078	1.46318	0.89286	0.59982	0.39919	0.32182	0.29057	0.24373	0.22776	0.18097	0.13408	0.14889	0.16298	0.16207	0.1613	0.08267	0.03562	0.03	0.02	0.01	0.005	0.005];

auxn=interp1(longdaexp_n*1e-9,n/10,lambda,'linear','extrap');
auxk=interp1(longdaexp_k*1e-9,k/10,lambda,'linear','extrap');

N=auxn+1i*auxk;

end

function N=PMMA(lambda)
%indice de refraccion del acrilico
% N. Sultanova, S. Kasarova and I. Nikolov. Dispersion properties of optical polymers, Acta Physica Polonica A 116, 585-587 (2009)
% (fit of the experimental data with the Sellmeier dispersion formula: Mikhail Polyanskiy)

lambda=lambda*1e6;
N=sqrt(1.1819*lambda.^2./(lambda.^2-0.011313)+1);

end

function N=Siamorfo(lambda)
%indice de refraccion del Silicio amorfo hidrogenado (queda por revisar en
%esta ecuacion en que unidades va lambda, probablemente en nanometros)
%buscar una referencia....
lambda=lambda*1e9;
n=3e5./lambda.^2+2.6;
alfa=10.^(1.5e6./lambda.^2-8);
k=lambda.*alfa/(4*pi);
N=n+1i*k;
end

function N=indiceSi5(lambda)
%funcion que devuelve el indice de refraccion real e imaginario del silicio
%puro. Interpola los valores de las tablas provistas por Sopra (ver
%http://refractiveindex.info/index.php?group=CRYSTALS&material=Si&option=so
%pra&wavelength=0.8266)

longda=[1.63138E-7	1.63569E-7	1.64002E-7	1.64437E-7	1.64874E-7	1.65314E-7	1.65756E-7	1.662E-7	1.66647E-7	1.67096E-7	1.67548E-7	1.68002E-7	1.68458E-7	1.68917E-7	1.69379E-7	1.69843E-7	1.70309E-7	1.70778E-7	1.7125E-7	1.71725E-7	1.72202E-7	1.72681E-7	1.73164E-7	1.73649E-7	1.74136E-7	1.74627E-7	1.7512E-7	1.75616E-7	1.76115E-7	1.76617E-7	1.77122E-7	1.77629E-7	1.7814E-7	1.78653E-7	1.79169E-7	1.79689E-7	1.80211E-7	1.80736E-7	1.81265E-7	1.81796E-7	1.82331E-7	1.82869E-7	1.8341E-7	1.83954E-7	1.84502E-7	1.85053E-7	1.85607E-7	1.86164E-7	1.86725E-7	1.87289E-7	1.87856E-7	1.88427E-7	1.89002E-7	1.8958E-7	1.90161E-7	1.90746E-7	1.91335E-7	1.91928E-7	1.92524E-7	1.93123E-7	1.93727E-7	1.94334E-7	1.94945E-7	1.9556E-7	1.96179E-7	1.96802E-7	1.97429E-7	1.98059E-7	1.98694E-7	1.99333E-7	1.99976E-7	2.00623E-7	2.01275E-7	2.0193E-7	2.0259E-7	2.03254E-7	2.03923E-7	2.04596E-7	2.05273E-7	2.05955E-7	2.06642E-7	2.07333E-7	2.08029E-7	2.08729E-7	2.09434E-7	2.10144E-7	2.10859E-7	2.11579E-7	2.12303E-7	2.13033E-7	2.13768E-7	2.14507E-7	2.15252E-7	2.16002E-7	2.16757E-7	2.17518E-7	2.18284E-7	2.19055E-7	2.19832E-7	2.20614E-7	2.21402E-7	2.22196E-7	2.22995E-7	2.238E-7	2.24611E-7	2.25428E-7	2.2625E-7	2.27079E-7	2.27914E-7	2.28755E-7	2.29602E-7	2.30456E-7	2.31316E-7	2.32182E-7	2.33055E-7	2.33934E-7	2.3482E-7	2.35713E-7	2.36613E-7	2.3752E-7	2.38433E-7	2.39354E-7	2.40281E-7	2.41216E-7	2.42159E-7	2.43108E-7	2.44065E-7	2.4503E-7	2.46002E-7	2.46982E-7	2.4797E-7	2.48966E-7	2.4997E-7	2.50982E-7	2.52002E-7	2.53031E-7	2.54068E-7	2.55114E-7	2.56168E-7	2.57231E-7	2.58302E-7	2.59383E-7	2.60473E-7	2.61572E-7	2.6268E-7	2.63798E-7	2.64926E-7	2.66063E-7	2.67209E-7	2.68366E-7	2.69533E-7	2.7071E-7	2.71897E-7	2.73095E-7	2.74303E-7	2.75523E-7	2.76753E-7	2.77994E-7	2.79246E-7	2.80509E-7	2.81785E-7	2.83071E-7	2.8437E-7	2.8568E-7	2.87003E-7	2.88338E-7	2.89685E-7	2.91045E-7	2.92418E-7	2.93804E-7	2.95203E-7	2.96615E-7	2.98041E-7	2.99481E-7	3.00935E-7	3.02403E-7	3.03885E-7	3.05382E-7	3.06894E-7	3.08421E-7	3.09963E-7	3.11521E-7	3.13094E-7	3.14683E-7	3.16289E-7	3.17911E-7	3.19549E-7	3.21205E-7	3.22878E-7	3.24569E-7	3.26277E-7	3.28003E-7	3.29748E-7	3.31511E-7	3.33293E-7	3.35095E-7	3.36916E-7	3.38757E-7	3.40619E-7	3.42501E-7	3.44403E-7	3.46327E-7	3.48273E-7	3.50241E-7	3.52231E-7	3.54243E-7	3.56279E-7	3.58339E-7	3.60422E-7	3.6253E-7	3.64662E-7	3.6682E-7	3.69004E-7	3.71213E-7	3.73449E-7	3.75713E-7	3.78004E-7	3.80323E-7	3.8267E-7	3.85047E-7	3.87454E-7	3.89891E-7	3.92358E-7	3.94857E-7	3.97388E-7	3.99952E-7	4.02549E-7	4.0518E-7	4.07846E-7	4.10547E-7	4.13284E-7	4.16058E-7	4.18869E-7	4.21718E-7	4.24607E-7	4.27535E-7	4.30504E-7	4.33515E-7	4.36568E-7	4.39664E-7	4.42804E-7	4.4599E-7	4.49222E-7	4.52501E-7	4.55828E-7	4.59204E-7	4.62631E-7	4.6611E-7	4.69641E-7	4.73226E-7	4.76866E-7	4.80563E-7	4.84317E-7	4.88131E-7	4.92005E-7	4.95941E-7	4.9994E-7	5.04005E-7	5.08136E-7	5.12335E-7	5.16605E-7	5.20946E-7	5.25361E-7	5.29851E-7	5.34419E-7	5.39066E-7	5.43795E-7	5.48607E-7	5.53505E-7	5.58492E-7	5.63569E-7	5.68739E-7	5.74005E-7	5.7937E-7	5.84836E-7	5.90406E-7	5.96083E-7	6.0187E-7	6.0777E-7	6.13788E-7	6.19926E-7	6.26188E-7	6.32577E-7	6.39099E-7	6.45756E-7	6.52554E-7	6.59496E-7	6.66587E-7	6.73833E-7	6.81237E-7	6.88807E-7	6.96546E-7	7.04461E-7	7.12559E-7	7.20844E-7	7.29325E-7	7.38007E-7	7.46899E-7	7.56007E-7	7.65341E-7	7.74907E-7	7.84716E-7	7.94777E-7	8.05099E-7	8.15692E-7	8.26568E-7	8.37738E-7	8.49214E-7	8.61008E-7	8.73135E-7	8.85608E-7	8.98443E-7	9.11656E-7	9.25263E-7	9.39282E-7	9.53732E-7	9.68634E-7	9.84009E-7	9.9988E-7	1.01627E-6	1.03321E-6	1.05072E-6	1.06884E-6	1.08759E-6	1.10701E-6	1.12714E-6	1.14801E-6	1.16967E-6	1.19217E-6	1.21554E-6	1.23985E-6	1.26515E-6	1.29151E-6	1.31899E-6	1.34767E-6	1.37761E-6	1.40892E-6	1.44169E-6	1.47601E-6	1.51201E-6	1.54981E-6	1.58955E-6	1.63138E-6	1.67548E-6	1.72202E-6	1.77122E-6	1.82331E-6	1.87856E-6	1.93727E-6	1.99976E-6	2.06642E-6	2.13768E-6	2.21402E-6	2.29602E-6	2.38433E-6	2.4797E-6 25e-6];
nSi=[0.52512	0.53294	0.54031	0.54727	0.55382	0.56	0.56583	0.57133	0.57652	0.58143	0.58608	0.59049	0.59469	0.5987	0.60254	0.60624	0.60982	0.61329	0.61669	0.62004	0.62336	0.62665	0.62996	0.63331	0.63675	0.64028	0.64389	0.64767	0.65164	0.65587	0.66045	0.66559	0.67109	0.67687	0.68287	0.68909	0.69549	0.70205	0.70874	0.71554	0.72238	0.72919	0.73602	0.7429	0.74982	0.75678	0.76373	0.77072	0.77778	0.78493	0.79217	0.79953	0.80702	0.81463	0.82238	0.83032	0.83851	0.84687	0.85537	0.86399	0.87271	0.88156	0.89048	0.89944	0.90839	0.91729	0.92603	0.93469	0.94333	0.95194	0.96052	0.969	0.97748	0.98636	0.99535	1.0035	1.00005	0.99525	0.99285	0.99656	1.01014	1.03522	1.04698	1.0655	1.0687	1.08037	1.09032	1.10141	1.10974	1.1188	1.13254	1.13911	1.15479	1.16458	1.17512	1.18589	1.19509	1.21069	1.22184	1.23471	1.24747	1.26525	1.28067	1.2984	1.31861	1.33997	1.36205	1.38926	1.41626	1.44469	1.47176	1.50138	1.52611	1.54821	1.56555	1.57855	1.58549	1.58987	1.5914	1.59189	1.5887	1.58614	1.58156	1.57829	1.57377	1.57107	1.56968	1.56882	1.56829	1.5689	1.57045	1.57495	1.57956	1.58427	1.59086	1.59704	1.60842	1.61756	1.62913	1.6432	1.65768	1.67327	1.69203	1.71318	1.73711	1.76411	1.7939	1.83134	1.8745	1.92745	1.98794	2.05927	2.14091	2.23435	2.33896	2.45173	2.57196	2.69982	2.83395	2.97449	3.12072	3.27712	3.44476	3.63446	3.85005	4.08699	4.31959	4.52564	4.6864	4.80472	4.88813	4.94119	4.97698	4.99945	5.01238	5.0196	5.02102	5.02014	5.0176	5.01459	5.01056	5.00866	5.0098	5.00902	5.01197	5.016	5.02156	5.02933	5.03963	5.05227	5.06495	5.07904	5.09542	5.11457	5.1345	5.15552	5.17913	5.20403	5.23097	5.26146	5.29565	5.33564	5.38341	5.44173	5.5149	5.60989	5.73358	5.89402	6.08828	6.30723	6.52147	6.69487	6.79634	6.8289	6.79788	6.70742	6.58566	6.45226	6.31582	6.18508	6.06231	5.94805	5.84163	5.74364	5.65376	5.56993	5.49236	5.42003	5.34916	5.28362	5.23397	5.17494	5.11914	5.06625	5.01599	4.9681	4.92238	4.87867	4.83684	4.79677	4.75835	4.72143	4.6859	4.65166	4.61864	4.58678	4.55608	4.52648	4.49797	4.47048	4.44393	4.41824	4.39332	4.36911	4.34556	4.32264	4.30034	4.27865	4.25754	4.237	4.21701	4.19754	4.17856	4.16005	4.14198	4.12435	4.10716	4.0904	4.07406	4.05814	4.04262	4.02748	4.0127	3.99829	3.98423	3.97051	3.9571	3.94399	3.93116	3.91857	3.90623	3.89307	3.88148	3.87026	3.85831	3.84675	3.83675	3.82573	3.8151	3.80537	3.79646	3.78732	3.77801	3.76841	3.76046	3.75218	3.74493	3.73629	3.72758	3.72091	3.71375	3.70525	3.69709	3.6877	3.68094	3.67186	3.6631	3.6548	3.64667	3.63873	3.6313	3.62484	3.6189	3.61328	3.60785	3.6025	3.597	3.59149	3.58602	3.58059	3.57517	3.5697	3.56424	3.55878	3.55334	3.54789	3.5424	3.53693	3.53151	3.52619	3.52099	3.51601	3.51117	3.50646	3.50186	3.49736	3.49292	3.48859	3.48438	3.48032	3.47641	3.47269	3.46914	3.46575	3.46254	3.4595	3.45664	3.45397	3.45141	3.44904	3.4471	3.44613	3.44548	3.44488	3.44406	3.44275 3.44275];
kSi=[2.16784	2.17386	2.18009	2.18652	2.19316	2.2	2.20703	2.21424	2.22164	2.22921	2.23696	2.24487	2.25294	2.26117	2.26955	2.27808	2.28675	2.29555	2.30448	2.31354	2.32272	2.33202	2.34143	2.35094	2.36055	2.37025	2.38005	2.38992	2.39986	2.40985	2.41987	2.42986	2.43986	2.44993	2.46006	2.47026	2.48053	2.4909	2.50137	2.51198	2.52274	2.53371	2.54485	2.55613	2.56755	2.57912	2.59084	2.60269	2.61465	2.6267	2.63882	2.65099	2.66321	2.67549	2.68781	2.7002	2.71271	2.72524	2.73774	2.75016	2.76253	2.77491	2.78712	2.79907	2.81064	2.82163	2.83164	2.84106	2.85011	2.85885	2.86723	2.87506	2.88277	2.89108	2.89978	2.90763	2.90534	2.90138	2.89873	2.90035	2.9092	2.92833	2.94214	2.93903	2.96087	2.97983	2.98963	3.00348	3.01565	3.02508	3.04478	3.06055	3.07263	3.08631	3.1027	3.11962	3.13499	3.15008	3.16915	3.19007	3.20674	3.22793	3.24533	3.26703	3.28519	3.30158	3.31862	3.33442	3.34957	3.35933	3.36639	3.36856	3.36746	3.36375	3.35851	3.35286	3.34632	3.34426	3.34369	3.34714	3.35363	3.36306	3.37556	3.38949	3.40823	3.42925	3.45084	3.47687	3.50371	3.53349	3.56474	3.59785	3.63254	3.67001	3.70933	3.74919	3.78895	3.83512	3.88012	3.92858	3.97909	4.03148	4.08838	4.14868	4.21115	4.2781	4.35062	4.42619	4.5063	4.59051	4.67808	4.76457	4.84967	4.93322	5.01133	5.08196	5.14776	5.20607	5.25773	5.30401	5.34382	5.38063	5.41369	5.43532	5.43846	5.39494	5.29987	5.15633	4.98743	4.81095	4.63919	4.47989	4.33429	4.20392	4.08585	3.97912	3.8845	3.79821	3.72005	3.64962	3.58645	3.52886	3.4767	3.42926	3.3855	3.34602	3.30958	3.27515	3.2418	3.21103	3.18167	3.15405	3.12746	3.10259	3.07897	3.05762	3.03873	3.02096	3.00651	2.99489	2.98661	2.98257	2.98356	2.98909	2.99908	3.01434	3.02674	3.02214	2.98097	2.88014	2.70489	2.45554	2.16783	1.86846	1.57646	1.31987	1.10942	0.94458	0.81485	0.71379	0.62979	0.56103	0.50462	0.45556	0.4164	0.38795	0.35583	0.32828	0.31238	0.29126	0.27769	0.26065	0.24528	0.23142	0.21895	0.20773	0.19757	0.1882	0.17937	0.17088	0.16275	0.15511	0.14822	0.14221	0.13705	0.1324	0.12774	0.12255	0.11654	0.10976	0.1026	0.09561	0.08926	0.08381	0.07928	0.07549	0.07217	0.06903	0.06588	0.06261	0.05924	0.05588	0.05268	0.04976	0.04715	0.04478	0.04255	0.04038	0.03826	0.03625	0.03441	0.03278	0.03136	0.03013	0.02901	0.02795	0.02686	0.02574	0.02457	0.02336	0.02213	0.02147	0.01954	0.01771	0.01678	0.01636	0.01572	0.01495	0.01404	0.01347	0.01304	0.01251	0.01191	0.01137	0.01086	0.01045	0.00994	0.00934	0.00877	0.00823	0.00772	0.00719	0.00663	0.00612	0.00561	0.00529	0.00498	0.00469	0.00443	0.0042	0.00399	0.00378	0.00358	0.0034	0.00322	0.00303	0.00283	0.00262	0.00242	0.00224	0.0021	0.00204	0.00202	0.00201	0.00198	0.00193	0.00179	0.00161	0.00142	0.00121	9.99E-4	7.78E-4	5.62E-4	3.54E-4	1.69E-4	3.1E-5	1E-6	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0  0 ];

aux=interp1(longda,nSi,lambda,'spline','extrap');
aux2=interp1(longda,kSi,lambda,'spline','extrap');
N=aux+1i*aux2;
end

function N=indiceSiO2(lambda)
%funcion que devuelve el indice de refraccion real e imaginario del cuarzo fundido.
%Interpola los valores de las tablas provistas por Sopra (ver
%http://refractiveindex.info/?group=GLASSES&material=F_SILICA

lambda=lambda*1e6; %longitud de onda en micrones
N = sqrt( 1 + 0.6961663*lambda.^2./(lambda.^2-power(0.0684043,2)) + 0.4079426*lambda.^2./(lambda.^2-power(0.1162414,2)) + 0.8974794*lambda.^2./(lambda.^2-power(9.896161,2)));

end

function N=indicePET(lambda)
%funcion que devuelve el indice de refraccion real e imaginario del PET (polietilentereftalato). Interpola los valores de las tablas provistas por Sopra (ver
%http://www.filmetrics.com/refractive-index-database/PET/Estar-Melinex-Mylar
% ver "Polarimetric characterization of optically anisotropic flexible substrates Thin Solid Films 516 (2008) 14141418"

longda=1e-9*[400 400.47	401.23	401.98	402.73	403.48	404.24	404.99	405.74	406.5	407.25	408	408.75	409.5	410.26	411.01	411.76	412.51	413.26	414.02	414.77	415.52	416.27	417.02	417.77	418.52	419.27	420.02	420.78	421.53	422.28	423.03	423.78	424.53	425.28	426.03	426.78	427.53	428.28	429.03	429.78	430.53	431.28	432.03	432.78	433.53	434.28	435.02	435.77	436.52	437.27	438.02	438.77	439.52	440.27	441.01	441.76	442.51	443.26	444.01	444.76	445.5	446.25	447	447.75	448.49	449.24	449.99	450.74	451.48	452.23	452.98	453.72	454.47	455.22	455.96	456.71	457.46	458.2	458.95	459.7	460.44	461.19	461.94	462.68	463.43	464.17	464.92	465.66	466.41	467.15	467.9	468.65	469.39	470.14	470.88	471.63	472.37	473.11	473.86	474.6	475.35	476.09	476.84	477.58	478.33	479.07	479.81	480.56	481.3	482.04	482.79	483.53	484.27	485.02	485.76	486.5	487.25	487.99	488.73	489.48	490.22	490.96	491.7	492.45	493.19	493.93	494.67	495.41	496.16	496.9	497.64	498.38	499.12	499.86	500.61	501.35	502.09	502.83	503.57	504.31	505.05	505.79	506.53	507.28	508.02	508.76	509.5	510.24	510.98	511.72	512.46	513.2	513.94	514.68	515.42	516.16	516.9	517.64	518.37	519.11	519.85	520.59	521.33	522.07	522.81	523.55	524.29	525.02	525.76	526.5	527.24	527.98	528.72	529.45	530.19	530.93	531.67	532.41	533.14	533.88	534.62	535.36	536.09	536.83	537.57	538.3	539.04	539.78	540.51	541.25	541.99	542.72	543.46	544.2	544.93	545.67	546.4	547.14	547.88	548.61	549.35	550.08	550.82	551.55	552.29	553.02	553.76	554.49	555.23	555.96	556.7	557.43	558.17	558.9	559.64	560.37	561.11	561.84	562.57	563.31	564.04	564.78	565.51	566.24	566.98	567.71	568.44	569.18	569.91	570.64	571.38	572.11	572.84	573.57	574.31	575.04	575.77	576.51	577.24	577.97	578.7	579.43	580.17	580.9	581.63	582.36	583.09	583.82	584.56	585.29	586.02	586.75	587.48	588.21	588.94	589.67	590.4	591.14	591.87	592.6	593.33	594.06	594.79	595.52	596.25	596.98	597.71	598.44	599.17	599.9	600.63	601.36	602.08	602.81	603.54	604.27	605	605.73	606.46	607.19	607.92	608.64	609.37	610.1	610.83	611.56	612.29	613.01	613.74	614.47	615.2	615.93	616.65	617.38	618.11	618.84	619.56	620.29	621.02	621.74	622.47	623.2	623.92	624.65	625.38	626.1	626.83	627.56	628.28	629.01	629.73	630.46	631.19	631.91	632.64	633.36	634.09	634.81	635.54	636.27	636.99	637.72	638.44	639.17	639.89	640.61	641.34	642.06	642.79	643.51	644.24	644.96	645.69	646.41	647.13	647.86	648.58	649.3	650.03	650.75	651.47	652.2	652.92	653.64	654.37	655.09	655.81	656.54	657.26	657.98	658.7	659.43	660.15	660.87	661.59	662.32	663.04	663.76	664.48	665.2	665.92	666.65	667.37	668.09	668.81	669.53	670.25	670.97	671.69	672.41	673.14	673.86	674.58	675.3	676.02	676.74	677.46	678.18	678.9	679.62	680.34	681.06	681.78	682.5	683.22	683.94	684.66	685.38	686.09	686.81	687.53	688.25	688.97	689.69	690.41	691.13	691.84	692.56	693.28	694	694.72	695.44	696.15	696.87	697.59	698.31	699.02	699.74	700];
n=[1.69295	1.6927	1.6922	1.6918	1.6914	1.691	1.6906	1.6902	1.6898	1.6894	1.689	1.6886	1.6882	1.6878	1.6874	1.687	1.6866	1.6863	1.6859	1.6855	1.6852	1.6848	1.6844	1.6841	1.6837	1.6834	1.683	1.6827	1.6823	1.682	1.6816	1.6813	1.681	1.6806	1.6803	1.68	1.6797	1.6793	1.679	1.6787	1.6784	1.6781	1.6778	1.6775	1.6771	1.6768	1.6765	1.6762	1.6759	1.6756	1.6754	1.6751	1.6748	1.6745	1.6742	1.6739	1.6736	1.6734	1.6731	1.6728	1.6725	1.6723	1.672	1.6717	1.6714	1.6712	1.6709	1.6707	1.6704	1.6701	1.6699	1.6696	1.6694	1.6691	1.6689	1.6686	1.6684	1.6681	1.6679	1.6676	1.6674	1.6672	1.6669	1.6667	1.6665	1.6662	1.666	1.6658	1.6655	1.6653	1.6651	1.6649	1.6646	1.6644	1.6642	1.664	1.6638	1.6635	1.6633	1.6631	1.6629	1.6627	1.6625	1.6623	1.6621	1.6619	1.6617	1.6615	1.6613	1.6611	1.6609	1.6607	1.6605	1.6603	1.6601	1.6599	1.6597	1.6595	1.6593	1.6591	1.6589	1.6587	1.6585	1.6584	1.6582	1.658	1.6578	1.6576	1.6575	1.6573	1.6571	1.6569	1.6567	1.6566	1.6564	1.6562	1.656	1.6559	1.6557	1.6555	1.6554	1.6552	1.655	1.6549	1.6547	1.6545	1.6544	1.6542	1.6541	1.6539	1.6537	1.6536	1.6534	1.6533	1.6531	1.6529	1.6528	1.6526	1.6525	1.6523	1.6522	1.652	1.6519	1.6517	1.6516	1.6514	1.6513	1.6512	1.651	1.6509	1.6507	1.6506	1.6504	1.6503	1.6502	1.65	1.6499	1.6497	1.6496	1.6495	1.6493	1.6492	1.6491	1.6489	1.6488	1.6487	1.6485	1.6484	1.6483	1.6481	1.648	1.6479	1.6477	1.6476	1.6475	1.6474	1.6472	1.6471	1.647	1.6469	1.6467	1.6466	1.6465	1.6464	1.6463	1.6461	1.646	1.6459	1.6458	1.6457	1.6455	1.6454	1.6453	1.6452	1.6451	1.645	1.6448	1.6447	1.6446	1.6445	1.6444	1.6443	1.6442	1.6441	1.6439	1.6438	1.6437	1.6436	1.6435	1.6434	1.6433	1.6432	1.6431	1.643	1.6429	1.6428	1.6427	1.6426	1.6425	1.6424	1.6423	1.6422	1.6421	1.642	1.6419	1.6418	1.6417	1.6416	1.6415	1.6414	1.6413	1.6412	1.6411	1.641	1.6409	1.6408	1.6407	1.6406	1.6405	1.6404	1.6403	1.6402	1.6401	1.64	1.6399	1.6398	1.6398	1.6397	1.6396	1.6395	1.6394	1.6393	1.6392	1.6391	1.639	1.639	1.6389	1.6388	1.6387	1.6386	1.6385	1.6384	1.6383	1.6383	1.6382	1.6381	1.638	1.6379	1.6378	1.6378	1.6377	1.6376	1.6375	1.6374	1.6374	1.6373	1.6372	1.6371	1.637	1.637	1.6369	1.6368	1.6367	1.6366	1.6366	1.6365	1.6364	1.6363	1.6363	1.6362	1.6361	1.636	1.636	1.6359	1.6358	1.6357	1.6357	1.6356	1.6355	1.6354	1.6354	1.6353	1.6352	1.6351	1.6351	1.635	1.6349	1.6349	1.6348	1.6347	1.6347	1.6346	1.6345	1.6344	1.6344	1.6343	1.6342	1.6342	1.6341	1.634	1.634	1.6339	1.6338	1.6338	1.6337	1.6336	1.6336	1.6335	1.6334	1.6334	1.6333	1.6332	1.6332	1.6331	1.6331	1.633	1.6329	1.6329	1.6328	1.6327	1.6327	1.6326	1.6326	1.6325	1.6324	1.6324	1.6323	1.6322	1.6322	1.6321	1.6321	1.632	1.6319	1.6319	1.6318	1.6318	1.6317	1.6317	1.6316	1.6315	1.6315	1.6314	1.6314	1.6313	1.6313	1.6312	1.6311	1.6311	1.631	1.631	1.6309	1.6309	1.6308	1.6308	1.6307	1.6306	1.6306	1.6305	1.6305	1.6304	1.6304	1.6303	1.6303	1.6302	1.6302	1.6301	1.6301	1.63	1.6299	1.6299];

aux=interp1(longda,n,lambda,'cubic','extrap');
N=aux;
end

function N=aluminio(lambda)
%indice de refraccion del aluminio obtenido de
%http://www-swiss.ai.mit.edu/~jaffer/FreeSnell/nk.html

longdaexp=[3.1e-5	2.48e-5	2.06667e-5	1.77143e-5	1.55e-5	1.37778e-5	1.24e-5	9.92e-6	8.26667e-6	7.08571e-6	6.2e-6	4.96e-6	4.13333e-6	3.54286e-6	3.1e-6	2.48e-6	2.06667e-6	1.77143e-6	1.55e-6	1.37778e-6	1.24e-6	1.12727e-6	1.03333e-6	9.53846e-7	8.85714e-7	8.26667e-7	7.75e-7	7.29412e-7	6.88889e-7	6.52632e-7	6.2e-7	5.63636e-7	5.16667e-7	4.76923e-7	4.42857e-7	4.13333e-7	3.875e-7	3.64706e-7	3.44444e-7	3.26316e-7	3.1e-7	2.95238e-7	2.81818e-7	2.69565e-7	2.58333e-7	2.48e-7	2.06667e-7	1.90769e-7	1.77143e-7	1.65333e-7	1.55e-7	1.45882e-7	1.37778e-7	1.30526e-7	1.24e-7	1.18095e-7	1.12727e-7	1.07826e-7	1.03333e-7	9.92e-8	9.53846e-8	9.18519e-8	8.85714e-8	8.73239e-8	8.61111e-8	8.49315e-8	8.37838e-8	8.26667e-8	8.15789e-8	8.05195e-8	7.94872e-8	7.8481e-8	7.75e-8	7.65432e-8	7.56098e-8	7.40299e-8	7.29412e-8	7.18841e-8	7.08571e-8	6.98592e-8	6.88889e-8	6.7027e-8	6.52632e-8	6.35897e-8	6.2e-8	6.04878e-8	5.90476e-8	5.76744e-8	5.63636e-8	5.51111e-8	5.3913e-8	5.2766e-8	5.16667e-8	5.06122e-8	4.96e-8	4.86275e-8	4.76923e-8	4.59259e-8	4.42857e-8	4.27586e-8	4.13333e-8	3.54286e-8	3.1e-8	2.75556e-8	2.48e-8	2.25455e-8	2.06667e-8	1.90769e-8	1.77143e-8	1.71034e-8	1.65333e-8	1.6e-8	1.55e-8	1.45882e-8	1.37778e-8	1.30526e-8	1.24e-8	1.12727e-8	1.03333e-8	9.53846e-9	8.85714e-9	8.26667e-9	7.75e-9	7.29412e-9	6.88889e-9	6.52632e-9	6.2e-9	5.63636e-9	5.16667e-9	4.76923e-9	4.42857e-9	4.13333e-9];
n=[98.595	74.997	62.852	53.79	45.784	39.651	34.464	24.965	18.572	14.274	11.733	8.586	6.759	5.438	4.454	3.072	2.273	1.77	1.444	1.264	1.212	1.201	1.26	1.468	2.237	2.745	2.625	2.143	1.741	1.488	1.304	1.018	0.826	0.695	0.598	0.523	0.46	0.407	0.363	0.326	0.294	0.267	0.244	0.223	0.205	0.19	0.13	0.11	0.095	0.082	0.072	0.063	0.056	0.049	0.044	0.04	0.036	0.033	0.033	0.034	0.038	0.041	0.048	0.053	0.058	0.067	0.086	0.125	0.178	0.234	0.28	0.318	0.351	0.38	0.407	0.448	0.474	0.498	0.52	0.54	0.558	0.591	0.62	0.646	0.668	0.689	0.707	0.724	0.739	0.753	0.766	0.778	0.789	0.799	0.809	0.817	0.826	0.84	0.854	0.865	0.876	0.915	0.94	0.957	0.969	0.979	0.987	0.995	1.006	1.025	1.011	1.008	1.007	1.007	1.005	0.999	0.991	0.994	0.991	0.987	0.989	0.99	0.989	0.989	0.99	0.99	0.991	0.992	0.993	0.993	0.994	0.995];
k=[203.701	172.199	150.799	135.5	123.734	114.102	105.6	89.25	76.96	66.93	59.37	48.235	40.96	35.599	31.485	25.581	21.403	18.328	15.955	14.021	12.464	11.181	10.01	8.949	8.212	8.309	8.597	8.573	8.205	7.821	7.479	6.846	6.283	5.8	5.385	5.024	4.708	4.426	4.174	3.946	3.74	3.552	3.38	3.222	3.076	2.942	2.391	2.173	1.983	1.814	1.663	1.527	1.402	1.286	1.178	1.076	0.979	0.883	0.791	0.7	0.609	0.517	0.417	0.373	0.327	0.273	0.211	0.153	0.108	0.184	0.073	0.065	0.06	0.055	0.05	0.045	0.042	0.04	0.038	0.036	0.035	0.032	0.03	0.028	0.027	0.025	0.024	0.023	0.022	0.021	0.021	0.02	0.019	0.018	0.018	0.017	0.016	0.015	0.014	0.014	0.013	0.01	0.008	0.007	0.006	0.005	0.004	0.004	0.004	0.004	0.024	0.025	0.024	0.028	0.031	0.036	0.03	0.025	0.024	0.021	0.016	0.015	0.014	0.011	0.01	0.009	0.007	0.006	0.005	0.004	0.003	0.002];


aux=interp1(longdaexp,n,lambda);
aux2=interp1(longdaexp,k,lambda);

N=aux+1i*aux2;
end

function N=aire(lambda)
%funcion que devuelve el indice de refraccion real e imaginario del aire
N=lambda*0+1.00029;
end

function N=BK7(lambda)
%indice de refraccion del vidrio obtenido de
%http://www-swiss.ai.mit.edu/~jaffer/FreeSnell/nk.html

longdaexp=[1.90746E-7	1.93727E-7	1.96802E-7	1.99976E-7	2.03254E-7	2.06642E-7	2.10144E-7	2.13768E-7	2.17518E-7	2.21402E-7	2.25428E-7	2.29602E-7	2.33934E-7	2.38433E-7	2.43108E-7	2.4797E-7	2.53031E-7	2.58302E-7	2.63798E-7	2.69533E-7	2.75523E-7	2.81784E-7	2.88338E-7	2.95203E-7	3.02403E-7	3.09963E-7	3.17911E-7	3.26277E-7	3.35095E-7	3.44403E-7	3.54243E-7	3.64662E-7	3.75713E-7	3.87454E-7	3.99952E-7	4.13284E-7	4.27535E-7	4.42804E-7	4.59204E-7	4.76866E-7	4.95941E-7	5.16605E-7	5.39066E-7	5.63569E-7	5.90406E-7	6.19926E-7	6.52554E-7	6.88807E-7	7.29325E-7	7.74907E-7	8.26568E-7	8.85608E-7	9.53732E-7	1.03321E-6	1.12714E-6	1.23985E-6];
n=[1.68543	1.67459	1.6645	1.65509	1.64632	1.63816	1.63056	1.62349	1.61689	1.61074	1.60499	1.59961	1.59457	1.58983	1.58538	1.58119	1.57722	1.57348	1.56993	1.56657	1.56337	1.56033	1.55743	1.55466	1.55201	1.54948	1.54705	1.54473	1.5425	1.54036	1.5383	1.53633	1.53443	1.53261	1.53085	1.52917	1.52755	1.526	1.5245	1.52307	1.52169	1.52036	1.51909	1.51786	1.51668	1.51554	1.51444	1.51337	1.51232	1.51129	1.51027	1.50924	1.50817	1.50705	1.50583	1.50444];

N=interp1(longdaexp,n,lambda);

end

function N=SnO2(lambda)
%indice de refraccion del SnO2 obtenido de
%Mater. Res. Soc. Symp. Proc., 426, (1996) 449
%se tomo el indice de refracion del oxido dopado con Fluor. Hay que tener
%mucho cuidado porque el indice del oxido puro es bastante diferente!!!

longdaexp=[3.0822E-7	3.1735E-7	3.2648E-7	3.3562E-7	3.4475E-7	3.5388E-7	3.6301E-7	3.7215E-7	3.8128E-7	3.9041E-7	3.9954E-7	4.0868E-7	4.1781E-7	4.2694E-7	4.3607E-7	4.4521E-7	4.5434E-7	4.6347E-7	4.726E-7	4.8174E-7	4.9087E-7	5E-7	5.0913E-7	5.1826E-7	5.274E-7	5.3653E-7	5.4566E-7	5.5479E-7	5.6393E-7	5.7306E-7	5.8219E-7	5.9132E-7	6.0046E-7	6.0959E-7	6.1872E-7	6.2785E-7	6.3699E-7	6.4612E-7	6.5525E-7	6.6438E-7	6.7352E-7	6.8265E-7	6.9178E-7	7.0091E-7	7.1005E-7	7.1918E-7	7.2831E-7	7.3744E-7	7.4658E-7	7.5571E-7	7.6484E-7	7.7397E-7	7.8311E-7	7.9224E-7	8.0137E-7	8.105E-7	8.1963E-7	8.2877E-7	8.379E-7	8.4703E-7	8.5616E-7	8.653E-7	8.7443E-7	8.8356E-7	8.9269E-7	9.0183E-7	9.1096E-7	9.2009E-7	9.2922E-7	9.3836E-7	9.4749E-7	9.5662E-7	9.6575E-7	9.7489E-7	9.8402E-7	9.9315E-7	1.0023E-6	1.0114E-6	1.0205E-6	1.0297E-6	1.0388E-6	1.0479E-6	1.0571E-6	1.0662E-6	1.0753E-6	1.0845E-6	1.0936E-6	1.1027E-6	1.1119E-6	1.121E-6	1.1301E-6	1.1393E-6	1.1484E-6	1.1575E-6	1.1667E-6	1.1758E-6	1.1849E-6	1.1941E-6	1.2032E-6	1.2123E-6	1.2215E-6	1.2306E-6	1.2397E-6	1.2489E-6	1.258E-6	1.2671E-6	1.2763E-6	1.2854E-6	1.2945E-6	1.3037E-6	1.3128E-6	1.3219E-6	1.3311E-6	1.3402E-6	1.3493E-6	1.3584E-6	1.3676E-6	1.3767E-6	1.3858E-6	1.395E-6	1.4041E-6	1.4132E-6	1.4224E-6	1.4315E-6	1.4406E-6	1.4498E-6	1.4589E-6	1.468E-6	1.4772E-6	1.4863E-6	1.4954E-6	1.5046E-6	1.5137E-6	1.5228E-6	1.532E-6	1.5411E-6	1.5502E-6	1.5594E-6	1.5685E-6	1.5776E-6	1.5868E-6	1.5959E-6	1.605E-6	1.6142E-6	1.6233E-6	1.6324E-6	1.6416E-6	1.6507E-6	1.6598E-6	1.6689E-6	1.6781E-6	1.6872E-6	1.6963E-6	1.7055E-6	1.7146E-6	1.7237E-6	1.7329E-6	1.742E-6	1.7511E-6	1.7603E-6	1.7694E-6	1.7785E-6	1.7877E-6	1.7968E-6	1.8059E-6	1.8151E-6	1.8242E-6	1.8333E-6	1.8425E-6	1.8516E-6	1.8607E-6	1.8699E-6	1.879E-6	1.8881E-6	1.8973E-6	1.9064E-6	1.9155E-6	1.9247E-6	1.9338E-6	1.9429E-6	1.9521E-6	1.9612E-6	1.9703E-6	1.9795E-6	1.9886E-6	1.9977E-6	2.0068E-6	2.016E-6	2.0251E-6	2.0342E-6	2.0434E-6	2.0525E-6	2.0616E-6	2.0708E-6	2.0799E-6	2.089E-6	2.0982E-6	2.1073E-6	2.1164E-6	2.1256E-6	2.1347E-6	2.1438E-6	2.153E-6	2.1621E-6	2.1712E-6	2.1804E-6	2.1895E-6	2.1986E-6	2.2078E-6	2.2169E-6	2.226E-6	2.2352E-6	2.2443E-6	2.2534E-6	2.2626E-6	2.2717E-6	2.2808E-6	2.29E-6	2.2991E-6	2.3082E-6	2.3174E-6	2.3265E-6	2.3356E-6	2.3447E-6	2.3539E-6	2.363E-6	2.3721E-6	2.3813E-6	2.3904E-6	2.3995E-6	2.4087E-6	2.4178E-6	2.4269E-6	2.4361E-6	2.4452E-6	2.4543E-6	2.4635E-6	2.4726E-6	2.4817E-6	2.4909E-6];
n=[2.2353	2.2206	2.2059	2.18238	2.15884	2.14118	2.12352	2.11174	2.09702	2.08526	2.0735	2.06762	2.05586	2.04704	2.04116	2.03234	2.02058	2.0147	2.00882	2	1.99412	1.98824	1.98236	1.97648	1.9706	1.96472	1.95884	1.95296	1.94708	1.94414	1.93826	1.93238	1.9265	1.92062	1.91474	1.90886	1.90592	1.90004	1.89416	1.88828	1.88238	1.87648	1.87058	1.86468	1.85878	1.85584	1.84996	1.84408	1.8382	1.83232	1.8235	1.81762	1.81174	1.80586	1.79998	1.7941	1.78528	1.7794	1.77352	1.76764	1.76176	1.75588	1.75	1.74412	1.7353	1.72648	1.7206	1.71472	1.7059	1.69708	1.6912	1.68532	1.6765	1.66768	1.6618	1.65592	1.64414	1.6353	1.6294	1.62056	1.60878	1.6029	1.59702	1.58526	1.57644	1.57056	1.56174	1.54998	1.5441	1.53528	1.52352	1.5147	1.50588	1.49412	1.4853	1.47648	1.46472	1.4559	1.44414	1.43238	1.42356	1.4118	1.40004	1.38826	1.37648	1.3647	1.35292	1.34114	1.32938	1.31762	1.29998	1.28822	1.27646	1.2647	1.25294	1.24118	1.22942	1.21178	1.20002	1.18238	1.17062	1.15296	1.14118	1.12352	1.11174	1.09408	1.07938	1.06174	1.0441	1.02646	1.00882	0.99412	0.97647	0.95882	0.94118	0.92353	0.9	0.88235	0.8647	0.84412	0.82647	0.80882	0.78529	0.76765	0.74706	0.72353	0.70588	0.68823	0.66471	0.64706	0.62941	0.61177	0.59412	0.57941	0.56177	0.55	0.53235	0.52059	0.50882	0.49706	0.48529	0.47647	0.46471	0.45588	0.45	0.44118	0.4353	0.43235	0.42647	0.42059	0.41764	0.4147	0.40882	0.40588	0.40294	0.4	0.39706	0.39412	0.39118	0.38823	0.38529	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38235	0.38529	0.38823	0.39118	0.39412	0.39706	0.39706	0.39706	0.39706	0.39706	0.39706	0.4	0.40294	0.40588	0.40882	0.41176	0.41176	0.41176	0.41176	0.41176	0.41176	0.41176	0.4147	0.41764	0.42059	0.42353	0.42647	0.42647	0.42647	0.42647];
k=[0.01229	0.01229	0.01229	0.01229	0.01229	0.01229	0.01229	0.01229	0.01229	0.01229	0.01229	0.01229	0.01229	0.01229	0.01229	0.01229	0.01229	0.01229	0.01229	0.01229	0.01229	0.01262	0.01309	0.0131	0.01312	0.01315	0.01319	0.01325	0.01332	0.0134	0.0135	0.01361	0.01373	0.01339	0.01319	0.01312	0.01317	0.01334	0.01363	0.01402	0.01451	0.01509	0.01576	0.0153	0.01585	0.01645	0.0171	0.01779	0.01852	0.01928	0.02006	0.02085	0.02166	0.02199	0.02244	0.02301	0.02366	0.02441	0.02522	0.02609	0.02701	0.02796	0.02892	0.0299	0.03087	0.03134	0.03191	0.03255	0.03324	0.03397	0.03471	0.03546	0.03618	0.03638	0.03725	0.0383	0.03951	0.04086	0.04232	0.04389	0.04553	0.04676	0.04815	0.04969	0.05134	0.05309	0.05442	0.05591	0.05754	0.05927	0.06106	0.06289	0.06424	0.0657	0.06781	0.07009	0.07201	0.07416	0.07648	0.07894	0.08103	0.08331	0.08526	0.08743	0.08978	0.09225	0.09492	0.09787	0.10058	0.10359	0.10685	0.10983	0.11258	0.11566	0.11958	0.12335	0.12701	0.13111	0.13508	0.13898	0.14332	0.14707	0.15193	0.15639	0.16155	0.16636	0.1719	0.17819	0.18475	0.19112	0.1979	0.20506	0.21211	0.21975	0.22808	0.23704	0.24611	0.25595	0.26615	0.27723	0.28861	0.30142	0.31459	0.32874	0.3434	0.35904	0.37552	0.39281	0.41042	0.42892	0.44826	0.46825	0.48884	0.5101	0.53196	0.5542	0.57675	0.59964	0.62279	0.64655	0.67036	0.69363	0.7184	0.7424	0.76672	0.79124	0.81584	0.83944	0.86347	0.88672	0.91026	0.93349	0.95582	0.97878	1.00131	1.02295	1.04465	1.0658	1.08751	1.10872	1.12896	1.14982	1.16964	1.18998	1.20982	1.22978	1.2493	1.2683	1.28729	1.30584	1.32509	1.3434	1.36179	1.3797	1.39828	1.41591	1.43361	1.45086	1.46784	1.48512	1.50155	1.51878	1.5358	1.55258	1.56905	1.5864	1.60242	1.61925	1.63588	1.6523	1.66847	1.68497	1.70068	1.71726	1.73253	1.74863	1.76458	1.78036	1.79596	1.81135	1.81135	1.81135	1.81135	1.81135	1.81135	1.81135	1.81135	1.81135	1.81135	1.81135	1.81135	1.81135	1.81135	1.81135	1.81135	1.81135	1.81135	1.81135	1.81135	1.81135];

aux=interp1(longdaexp,n,lambda);
aux2=interp1(longdaexp,k,lambda);

N=aux+aux2*1i;
end

function N=TiO2film(lambda)
%indice de refraccion del Titanium(IV) oxide (Titanium dioxide, TiO2) thin film (thickness 200 nm) on glass substrate.
%S. Sarkar, V. Gupta, M. Kumar, J. Schubert, P.T. Probst, J. Joseph, T.A.F. König, Hybridized guided-mode resonances via colloidal plasmonic self-assembled grating, ACS Appl. Mater. Interfaces, 11, 13752-13760 (2019)
%ver https://refractiveindex.info/?shelf=main&book=TiO2&page=Sarkar

longdaexp=[0.3	0.301	0.302	0.303	0.304	0.305	0.306	0.307	0.308	0.309	0.31	0.311	0.312	0.313	0.314	0.315	0.316	0.317	0.318	0.319	0.32	0.321	0.322	0.323	0.324	0.325	0.326	0.327	0.328	0.329	0.33	0.331	0.332	0.333	0.334	0.335	0.336	0.337	0.338	0.339	0.34	0.341	0.342	0.343	0.344	0.345	0.346	0.347	0.348	0.349	0.35	0.351	0.352	0.353	0.354	0.355	0.356	0.357	0.358	0.359	0.36	0.361	0.362	0.363	0.364	0.365	0.366	0.367	0.368	0.369	0.37	0.371	0.372	0.373	0.374	0.375	0.376	0.377	0.378	0.379	0.38	0.381	0.382	0.383	0.384	0.385	0.386	0.387	0.388	0.389	0.39	0.391	0.392	0.393	0.394	0.395	0.396	0.397	0.398	0.399	0.4	0.401	0.402	0.403	0.404	0.405	0.406	0.407	0.408	0.409	0.41	0.411	0.412	0.413	0.414	0.415	0.416	0.417	0.418	0.419	0.42	0.421	0.422	0.423	0.424	0.425	0.426	0.427	0.428	0.429	0.43	0.431	0.432	0.433	0.434	0.435	0.436	0.437	0.438	0.439	0.44	0.441	0.442	0.443	0.444	0.445	0.446	0.447	0.448	0.449	0.45	0.451	0.452	0.453	0.454	0.455	0.456	0.457	0.458	0.459	0.46	0.461	0.462	0.463	0.464	0.465	0.466	0.467	0.468	0.469	0.47	0.471	0.472	0.473	0.474	0.475	0.476	0.477	0.478	0.479	0.48	0.481	0.482	0.483	0.484	0.485	0.486	0.487	0.488	0.489	0.49	0.491	0.492	0.493	0.494	0.495	0.496	0.497	0.498	0.499	0.5	0.501	0.502	0.503	0.504	0.505	0.506	0.507	0.508	0.509	0.51	0.511	0.512	0.513	0.514	0.515	0.516	0.517	0.518	0.519	0.52	0.521	0.522	0.523	0.524	0.525	0.526	0.527	0.528	0.529	0.53	0.531	0.532	0.533	0.534	0.535	0.536	0.537	0.538	0.539	0.54	0.541	0.542	0.543	0.544	0.545	0.546	0.547	0.548	0.549	0.55	0.551	0.552	0.553	0.554	0.555	0.556	0.557	0.558	0.559	0.56	0.561	0.562	0.563	0.564	0.565	0.566	0.567	0.568	0.569	0.57	0.571	0.572	0.573	0.574	0.575	0.576	0.577	0.578	0.579	0.58	0.581	0.582	0.583	0.584	0.585	0.586	0.587	0.588	0.589	0.59	0.591	0.592	0.593	0.594	0.595	0.596	0.597	0.598	0.599	0.6	0.601	0.602	0.603	0.604	0.605	0.606	0.607	0.608	0.609	0.61	0.611	0.612	0.613	0.614	0.615	0.616	0.617	0.618	0.619	0.62	0.621	0.622	0.623	0.624	0.625	0.626	0.627	0.628	0.629	0.63	0.631	0.632	0.633	0.634	0.635	0.636	0.637	0.638	0.639	0.64	0.641	0.642	0.643	0.644	0.645	0.646	0.647	0.648	0.649	0.65	0.651	0.652	0.653	0.654	0.655	0.656	0.657	0.658	0.659	0.66	0.661	0.662	0.663	0.664	0.665	0.666	0.667	0.668	0.669	0.67	0.671	0.672	0.673	0.674	0.675	0.676	0.677	0.678	0.679	0.68	0.681	0.682	0.683	0.684	0.685	0.686	0.687	0.688	0.689	0.69	0.691	0.692	0.693	0.694	0.695	0.696	0.697	0.698	0.699	0.7	0.701	0.702	0.703	0.704	0.705	0.706	0.707	0.708	0.709	0.71	0.711	0.712	0.713	0.714	0.715	0.716	0.717	0.718	0.719	0.72	0.721	0.722	0.723	0.724	0.725	0.726	0.727	0.728	0.729	0.73	0.731	0.732	0.733	0.734	0.735	0.736	0.737	0.738	0.739	0.74	0.741	0.742	0.743	0.744	0.745	0.746	0.747	0.748	0.749	0.75	0.751	0.752	0.753	0.754	0.755	0.756	0.757	0.758	0.759	0.76	0.761	0.762	0.763	0.764	0.765	0.766	0.767	0.768	0.769	0.77	0.771	0.772	0.773	0.774	0.775	0.776	0.777	0.778	0.779	0.78	0.781	0.782	0.783	0.784	0.785	0.786	0.787	0.788	0.789	0.79	0.791	0.792	0.793	0.794	0.795	0.796	0.797	0.798	0.799	0.8	0.801	0.802	0.803	0.804	0.805	0.806	0.807	0.808	0.809	0.81	0.811	0.812	0.813	0.814	0.815	0.816	0.817	0.818	0.819	0.82	0.821	0.822	0.823	0.824	0.825	0.826	0.827	0.828	0.829	0.83	0.831	0.832	0.833	0.834	0.835	0.836	0.837	0.838	0.839	0.84	0.841	0.842	0.843	0.844	0.845	0.846	0.847	0.848	0.849	0.85	0.851	0.852	0.853	0.854	0.855	0.856	0.857	0.858	0.859	0.86	0.861	0.862	0.863	0.864	0.865	0.866	0.867	0.868	0.869	0.87	0.871	0.872	0.873	0.874	0.875	0.876	0.877	0.878	0.879	0.88	0.881	0.882	0.883	0.884	0.885	0.886	0.887	0.888	0.889	0.89	0.891	0.892	0.893	0.894	0.895	0.896	0.897	0.898	0.899	0.9	0.901	0.902	0.903	0.904	0.905	0.906	0.907	0.908	0.909	0.91	0.911	0.912	0.913	0.914	0.915	0.916	0.917	0.918	0.919	0.92	0.921	0.922	0.923	0.924	0.925	0.926	0.927	0.928	0.929	0.93	0.931	0.932	0.933	0.934	0.935	0.936	0.937	0.938	0.939	0.94	0.941	0.942	0.943	0.944	0.945	0.946	0.947	0.948	0.949	0.95	0.951	0.952	0.953	0.954	0.955	0.956	0.957	0.958	0.959	0.96	0.961	0.962	0.963	0.964	0.965	0.966	0.967	0.968	0.969	0.97	0.971	0.972	0.973	0.974	0.975	0.976	0.977	0.978	0.979	0.98	0.981	0.982	0.983	0.984	0.985	0.986	0.987	0.988	0.989	0.99	0.991	0.992	0.993	0.994	0.995	0.996	0.997	0.998	0.999	1	1.0025	1.005	1.0075	1.01	1.0125	1.015	1.0175	1.02	1.0225	1.025	1.0275	1.03	1.0325	1.035	1.0375	1.04	1.0425	1.045	1.0475	1.05	1.0525	1.055	1.0575	1.06	1.0625	1.065	1.0675	1.07	1.0725	1.075	1.0775	1.08	1.0825	1.085	1.0875	1.09	1.0925	1.095	1.0975	1.1	1.1025	1.105	1.1075	1.11	1.1125	1.115	1.1175	1.12	1.1225	1.125	1.1275	1.13	1.1325	1.135	1.1375	1.14	1.1425	1.145	1.1475	1.15	1.1525	1.155	1.1575	1.16	1.1625	1.165	1.1675	1.17	1.1725	1.175	1.1775	1.18	1.1825	1.185	1.1875	1.19	1.1925	1.195	1.1975	1.2	1.2025	1.205	1.2075	1.21	1.2125	1.215	1.2175	1.22	1.2225	1.225	1.2275	1.23	1.2325	1.235	1.2375	1.24	1.2425	1.245	1.2475	1.25	1.2525	1.255	1.2575	1.26	1.2625	1.265	1.2675	1.27	1.2725	1.275	1.2775	1.28	1.2825	1.285	1.2875	1.29	1.2925	1.295	1.2975	1.3	1.3025	1.305	1.3075	1.31	1.3125	1.315	1.3175	1.32	1.3225	1.325	1.3275	1.33	1.3325	1.335	1.3375	1.34	1.3425	1.345	1.3475	1.35	1.3525	1.355	1.3575	1.36	1.3625	1.365	1.3675	1.37	1.3725	1.375	1.3775	1.38	1.3825	1.385	1.3875	1.39	1.3925	1.395	1.3975	1.4	1.4025	1.405	1.4075	1.41	1.4125	1.415	1.4175	1.42	1.4225	1.425	1.4275	1.43	1.4325	1.435	1.4375	1.44	1.4425	1.445	1.4475	1.45	1.4525	1.455	1.4575	1.46	1.4625	1.465	1.4675	1.47	1.4725	1.475	1.4775	1.48	1.4825	1.485	1.4875	1.49	1.4925	1.495	1.4975	1.5	1.5025	1.505	1.5075	1.51	1.5125	1.515	1.5175	1.52	1.5225	1.525	1.5275	1.53	1.5325	1.535	1.5375	1.54	1.5425	1.545	1.5475	1.55	1.5525	1.555	1.5575	1.56	1.5625	1.565	1.5675	1.57	1.5725	1.575	1.5775	1.58	1.5825	1.585	1.5875	1.59	1.5925	1.595	1.5975	1.6	1.6025	1.605	1.6075	1.61	1.6125	1.615	1.6175	1.62	1.6225	1.625	1.6275	1.63	1.6325	1.635	1.6375	1.64	1.6425	1.645	1.6475	1.65	1.6525	1.655	1.6575	1.66	1.6625	1.665	1.6675	1.67	1.6725	1.675	1.6775	1.68	1.6825	1.685	1.6875	1.699]*1e-6;
n=[2.80998	2.81342	2.81645	2.81907	2.82127	2.82307	2.82445	2.82542	2.82599	2.82615	2.82591	2.82527	2.82424	2.82282	2.82102	2.81885	2.8163	2.8134	2.81014	2.80655	2.80261	2.79836	2.79378	2.78891	2.78374	2.77829	2.77257	2.76658	2.76035	2.75389	2.74719	2.74029	2.73319	2.7259	2.71843	2.7108	2.70302	2.69511	2.68707	2.67891	2.67065	2.66231	2.65389	2.64541	2.63687	2.6283	2.6197	2.61108	2.60246	2.59386	2.58527	2.57672	2.56822	2.55977	2.5514	2.54312	2.53494	2.52688	2.51894	2.51116	2.50354	2.49611	2.4889	2.48194	2.4753	2.4691	2.46325	2.45764	2.45227	2.44708	2.44208	2.43723	2.43253	2.42798	2.42355	2.41924	2.41505	2.41096	2.40698	2.40309	2.3993	2.39559	2.39197	2.38843	2.38497	2.38157	2.37825	2.375	2.37181	2.36869	2.36563	2.36262	2.35968	2.35678	2.35394	2.35115	2.34842	2.34572	2.34308	2.34048	2.33793	2.33542	2.33295	2.33051	2.32812	2.32577	2.32345	2.32117	2.31893	2.31672	2.31454	2.3124	2.31028	2.3082	2.30615	2.30413	2.30213	2.30017	2.29823	2.29632	2.29444	2.29258	2.29074	2.28894	2.28715	2.28539	2.28365	2.28194	2.28024	2.27857	2.27692	2.27529	2.27368	2.27209	2.27052	2.26897	2.26744	2.26592	2.26443	2.26295	2.26149	2.26005	2.25862	2.25721	2.25582	2.25444	2.25308	2.25173	2.2504	2.24909	2.24778	2.24649	2.24522	2.24396	2.24272	2.24148	2.24026	2.23906	2.23786	2.23668	2.23551	2.23436	2.23321	2.23208	2.23096	2.22985	2.22875	2.22766	2.22659	2.22552	2.22447	2.22342	2.22239	2.22136	2.22035	2.21934	2.21835	2.21736	2.21639	2.21542	2.21446	2.21351	2.21257	2.21164	2.21072	2.2098	2.2089	2.208	2.20711	2.20623	2.20536	2.20449	2.20364	2.20279	2.20194	2.20111	2.20028	2.19946	2.19865	2.19784	2.19704	2.19625	2.19547	2.19469	2.19392	2.19315	2.19239	2.19164	2.19089	2.19015	2.18942	2.18869	2.18797	2.18725	2.18654	2.18584	2.18514	2.18444	2.18376	2.18308	2.1824	2.18173	2.18106	2.1804	2.17975	2.1791	2.17845	2.17781	2.17718	2.17655	2.17592	2.1753	2.17468	2.17407	2.17347	2.17286	2.17227	2.17168	2.17109	2.1705	2.16992	2.16935	2.16878	2.16821	2.16765	2.16709	2.16654	2.16599	2.16544	2.1649	2.16436	2.16382	2.16329	2.16277	2.16224	2.16172	2.16121	2.1607	2.16019	2.15968	2.15918	2.15868	2.15819	2.1577	2.15721	2.15672	2.15624	2.15577	2.15529	2.15482	2.15435	2.15389	2.15343	2.15297	2.15251	2.15206	2.15161	2.15116	2.15072	2.15028	2.14984	2.14941	2.14898	2.14855	2.14812	2.1477	2.14728	2.14686	2.14644	2.14603	2.14562	2.14521	2.14481	2.14441	2.14401	2.14361	2.14322	2.14282	2.14243	2.14205	2.14166	2.14128	2.1409	2.14052	2.14015	2.13977	2.1394	2.13903	2.13867	2.1383	2.13794	2.13758	2.13723	2.13687	2.13652	2.13617	2.13582	2.13547	2.13513	2.13478	2.13444	2.1341	2.13377	2.13343	2.1331	2.13277	2.13244	2.13212	2.13179	2.13147	2.13115	2.13083	2.13051	2.13019	2.12988	2.12957	2.12926	2.12895	2.12864	2.12834	2.12804	2.12773	2.12743	2.12714	2.12684	2.12655	2.12625	2.12596	2.12567	2.12538	2.1251	2.12481	2.12453	2.12425	2.12397	2.12369	2.12341	2.12313	2.12286	2.12259	2.12232	2.12205	2.12178	2.12151	2.12125	2.12098	2.12072	2.12046	2.1202	2.11994	2.11968	2.11943	2.11917	2.11892	2.11867	2.11842	2.11817	2.11792	2.11767	2.11743	2.11719	2.11694	2.1167	2.11646	2.11622	2.11598	2.11575	2.11551	2.11528	2.11505	2.11481	2.11458	2.11436	2.11413	2.1139	2.11367	2.11345	2.11322	2.113	2.11278	2.11256	2.11234	2.11212	2.11191	2.11169	2.11148	2.11126	2.11105	2.11084	2.11063	2.11042	2.11021	2.11	2.1098	2.10959	2.10939	2.10918	2.10898	2.10878	2.10858	2.10838	2.10818	2.10798	2.10778	2.10759	2.10739	2.1072	2.10701	2.10682	2.10662	2.10643	2.10624	2.10606	2.10587	2.10568	2.1055	2.10531	2.10513	2.10494	2.10476	2.10458	2.1044	2.10422	2.10404	2.10386	2.10368	2.10351	2.10333	2.10316	2.10298	2.10281	2.10264	2.10247	2.10229	2.10212	2.10195	2.10179	2.10162	2.10145	2.10128	2.10112	2.10095	2.10079	2.10063	2.10046	2.1003	2.10014	2.09998	2.09982	2.09966	2.0995	2.09934	2.09918	2.09903	2.09887	2.09872	2.09856	2.09841	2.09826	2.0981	2.09795	2.0978	2.09765	2.0975	2.09735	2.0972	2.09706	2.09691	2.09676	2.09662	2.09647	2.09633	2.09618	2.09604	2.0959	2.09575	2.09561	2.09547	2.09533	2.09519	2.09505	2.09491	2.09477	2.09464	2.0945	2.09436	2.09423	2.09409	2.09396	2.09382	2.09369	2.09356	2.09342	2.09329	2.09316	2.09303	2.0929	2.09277	2.09264	2.09251	2.09238	2.09226	2.09213	2.092	2.09188	2.09175	2.09163	2.0915	2.09138	2.09125	2.09113	2.09101	2.09088	2.09076	2.09064	2.09052	2.0904	2.09028	2.09016	2.09004	2.08993	2.08981	2.08969	2.08957	2.08946	2.08934	2.08922	2.08911	2.089	2.08888	2.08877	2.08865	2.08854	2.08843	2.08832	2.0882	2.08809	2.08798	2.08787	2.08776	2.08765	2.08754	2.08744	2.08733	2.08722	2.08711	2.08701	2.0869	2.08679	2.08669	2.08658	2.08648	2.08637	2.08627	2.08616	2.08606	2.08596	2.08585	2.08575	2.08565	2.08555	2.08545	2.08535	2.08525	2.08515	2.08505	2.08495	2.08485	2.08475	2.08465	2.08455	2.08446	2.08436	2.08426	2.08417	2.08407	2.08398	2.08388	2.08379	2.08369	2.0836	2.0835	2.08341	2.08332	2.08322	2.08313	2.08304	2.08295	2.08285	2.08276	2.08267	2.08258	2.08249	2.0824	2.08231	2.08222	2.08213	2.08205	2.08196	2.08187	2.08178	2.08169	2.08161	2.08152	2.08143	2.08135	2.08126	2.08118	2.08109	2.08101	2.08092	2.08084	2.08075	2.08067	2.08059	2.0805	2.08042	2.08034	2.08026	2.08017	2.08009	2.08001	2.07993	2.07985	2.07977	2.07969	2.07961	2.07953	2.07945	2.07937	2.07929	2.07921	2.07913	2.07905	2.07897	2.0789	2.07882	2.07874	2.07867	2.07859	2.07851	2.07844	2.07836	2.07829	2.07821	2.07814	2.07806	2.07799	2.07791	2.07784	2.07776	2.07769	2.07762	2.07754	2.07747	2.0774	2.07733	2.07725	2.07718	2.07711	2.07704	2.07697	2.07689	2.07682	2.07675	2.07668	2.07661	2.07654	2.07647	2.0764	2.07634	2.07627	2.0762	2.07613	2.07606	2.07599	2.07593	2.07586	2.07579	2.07572	2.07566	2.07549	2.07532	2.07516	2.075	2.07484	2.07468	2.07452	2.07436	2.07421	2.07405	2.0739	2.07374	2.07359	2.07344	2.07329	2.07315	2.073	2.07285	2.07271	2.07257	2.07242	2.07228	2.07214	2.072	2.07186	2.07173	2.07159	2.07146	2.07132	2.07119	2.07106	2.07093	2.0708	2.07067	2.07054	2.07041	2.07028	2.07016	2.07003	2.06991	2.06979	2.06967	2.06955	2.06942	2.0693	2.06919	2.06907	2.06895	2.06884	2.06872	2.06861	2.06849	2.06838	2.06827	2.06816	2.06805	2.06794	2.06783	2.06772	2.06761	2.06751	2.0674	2.0673	2.06719	2.06709	2.06698	2.06688	2.06678	2.06668	2.06658	2.06648	2.06638	2.06628	2.06618	2.06609	2.06599	2.06589	2.0658	2.0657	2.06561	2.06552	2.06542	2.06533	2.06524	2.06515	2.06506	2.06497	2.06488	2.06479	2.06471	2.06462	2.06453	2.06444	2.06436	2.06427	2.06419	2.0641	2.06402	2.06394	2.06386	2.06377	2.06369	2.06361	2.06353	2.06345	2.06337	2.06329	2.06321	2.06313	2.06306	2.06298	2.0629	2.06283	2.06275	2.06268	2.0626	2.06253	2.06245	2.06238	2.0623	2.06223	2.06216	2.06209	2.06202	2.06195	2.06188	2.06181	2.06174	2.06167	2.0616	2.06153	2.06146	2.06139	2.06133	2.06126	2.06119	2.06113	2.06106	2.061	2.06093	2.06087	2.0608	2.06074	2.06068	2.06061	2.06055	2.06049	2.06043	2.06036	2.0603	2.06024	2.06018	2.06012	2.06006	2.06	2.05994	2.05988	2.05982	2.05977	2.05971	2.05965	2.05959	2.05954	2.05948	2.05942	2.05937	2.05931	2.05925	2.0592	2.05915	2.05909	2.05904	2.05898	2.05893	2.05888	2.05882	2.05877	2.05872	2.05866	2.05861	2.05856	2.05851	2.05846	2.05841	2.05836	2.05831	2.05826	2.05821	2.05816	2.05811	2.05806	2.05801	2.05796	2.05791	2.05786	2.05782	2.05777	2.05772	2.05767	2.05763	2.05758	2.05754	2.05749	2.05744	2.0574	2.05735	2.05731	2.05726	2.05722	2.05717	2.05713	2.05708	2.05704	2.057	2.05695	2.05691	2.05687	2.05682	2.05678	2.05674	2.0567	2.05666	2.05661	2.05657	2.05653	2.05649	2.05645	2.05641	2.05637	2.05633	2.05629	2.05625	2.05621	2.05617	2.05613	2.05609	2.05605	2.05601	2.05597	2.05594	2.0559	2.05586	2.05582	2.05578	2.05575	2.05571	2.05567	2.05564	2.0556	2.05556	2.05553	2.05549	2.05545	2.05542	2.05538	2.05535	2.05531	2.05528	2.05524	2.05521	2.05517	2.05514	2.0551	2.05507	2.05503	2.055	2.05497	2.05493	2.0549	2.05486	2.05483	2.0548	2.05477	2.05473	2.0547	2.054670];
k=[0.59278	0.57775	0.56267	0.54755	0.53241	0.51728	0.50215	0.48705	0.472	0.457	0.44208	0.42724	0.41251	0.3979	0.38342	0.36907	0.35489	0.34088	0.32705	0.31342	0.29999	0.28678	0.2738	0.26106	0.24856	0.23632	0.22434	0.21264	0.20121	0.19007	0.17923	0.16868	0.15843	0.14849	0.13887	0.12956	0.12057	0.1119	0.10356	0.09554	0.08786	0.0805	0.07347	0.06677	0.0604	0.05436	0.04865	0.04327	0.03822	0.03349	0.02909	0.025	0.02124	0.0178	0.01467	0.01185	0.00935	0.00715	0.00525	0.00365	0.00234	0.00133	6.07E-4	1.65E-4	1E-6	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];

aux=interp1(longdaexp,n,lambda);
aux2=interp1(longdaexp,k,lambda);

N=aux+aux2*1i;
end


function N = MAPI_Forouhi(lambda)
% Función que devuelve el índice de refracción y el coeficiente de
% extinción de MAPI, encontrados con el MODELO DE FOROUHI-BLOOMER con 3 picos. Los
% parámetros fueron obtenidos, a partir de el mejor ajuste del n y k de
% Lin et al 2014 (función anterior a ésta)
% 
% En nm
E=1240e-9./lambda;

%Modelo de Forouhi-Bloomer con 3 picos 
Par1=[   1.9083514e-01   3.1130047e+00   2.4296540e+00   1.5458414e+00];
Par2=[   5.0683217e-02    4.9149622e+00   6.0505075e+00   2.2911040e+00];
Par3=[   3.2804407e+00    3.1376511e+00   7.7245215e+00   1.9263747e+00];

%pico principal (1)
A1=Par1(1);
B1=Par1(2);
C1=Par1(3);
Eg1=Par1(4);

k_FB1=(A1*(E-Eg1).^2)./(E.^2-B1*E+C1);
k_FB1=k_FB1.*(E>Eg1);

%pico de fonfo (3)

A3=Par3(1);
B3=Par3(2);
C3=Par3(3);
Eg3=Par3(4);

k_FB3=(A3*(E-Eg3).^2)./(E.^2-B3*E+C3);
k_FB3=k_FB3.*(E>Eg3);

%pico secundario (2)

A2=Par2(1);
B2=Par2(2);
C2=Par2(3);
Eg2=Par2(4);

k_FB2=(A2*(E-Eg2).^2)./(E.^2-B2*E+C2);
k_FB2=k_FB2.*(E>Eg2);

% el k total es la suma de los tres picos
k=k_FB1+k_FB3+k_FB2;

%calculo la parte real
n_inf=2.8;

%obtenemos n a partir de Kramers-Kronig 

%pico principal
Q1=1/2*sqrt(4*C1-B1^2);
B01=A1/Q1*(-B1^2/2+B1*Eg1-Eg1^2+C1);
C01=A1/Q1*((Eg1^2+C1)*B1/2-2*Eg1*C1);

n1=(B01*E+C01)./(E.^2-B1*E+C1);

% pico de fondo
Q3=1/2*sqrt(4*C3-B3^2);
B03=A3/Q3*(-B3^2/2+B3*Eg3-Eg3^2+C3);
C03=A3/Q3*((Eg3^2+C3)*B3/2-2*Eg3*C3);

n3=(B03*E+C03)./(E.^2-B3*E+C3);

% pico secundario

Q2=1/2*sqrt(4*C2-B2^2);
B02=A2/Q2*(-B2^2/2+B2*Eg2-Eg2^2+C2);
C02=A2/Q2*((Eg2^2+C2)*B2/2-2*Eg2*C2);

n2=(B02*E+C02)./(E.^2-B2*E+C2);

%le sumo la parte continua
n=n1+n2+n3+n_inf;

N=n+1i*k;

end



function N=PbI2(lambda)
%indice de refraccion del PbI2
%ver https://refractiveindex.info/?shelf=main&book=PbI2&page=Frisenda 

longdaexp=[0.42346	0.42383	0.4242	0.42457	0.42494	0.42531	0.42568	0.42605	0.42642	0.42679	0.42716	0.42753	0.4279	0.42827	0.42864	0.42901	0.42938	0.42975	0.43012	0.43049	0.43086	0.43123	0.4316	0.43197	0.43234	0.43271	0.43308	0.43345	0.43382	0.43419	0.43456	0.43493	0.4353	0.43567	0.43604	0.43641	0.43678	0.43715	0.43752	0.43789	0.43826	0.43863	0.439	0.43937	0.43974	0.44011	0.44048	0.44085	0.44122	0.44159	0.44196	0.44233	0.4427	0.44307	0.44344	0.44381	0.44418	0.44455	0.44492	0.44529	0.44566	0.44603	0.4464	0.44677	0.44714	0.44751	0.44788	0.44825	0.44862	0.44899	0.44936	0.44973	0.4501	0.45047	0.45084	0.45121	0.45158	0.45195	0.45232	0.4527	0.45307	0.45344	0.45381	0.45418	0.45455	0.45492	0.45529	0.45566	0.45603	0.4564	0.45677	0.45714	0.45751	0.45788	0.45825	0.45862	0.45899	0.45936	0.45973	0.4601	0.46047	0.46084	0.46121	0.46158	0.46195	0.46232	0.46269	0.46306	0.46343	0.4638	0.46417	0.46454	0.46491	0.46528	0.46565	0.46602	0.46639	0.46676	0.46713	0.4675	0.46787	0.46824	0.46861	0.46898	0.46935	0.46972	0.47009	0.47046	0.47083	0.4712	0.47157	0.47194	0.47231	0.47268	0.47305	0.47342	0.47379	0.47416	0.47453	0.4749	0.47527	0.47564	0.47601	0.47638	0.47675	0.47712	0.47749	0.47786	0.47823	0.4786	0.47897	0.47934	0.47971	0.48008	0.48045	0.48082	0.48119	0.48156	0.48193	0.4823	0.48267	0.48304	0.48341	0.48378	0.48415	0.48452	0.48489	0.48526	0.48563	0.486	0.48637	0.48674	0.48711	0.48748	0.48785	0.48822	0.48859	0.48896	0.48933	0.4897	0.49007	0.49044	0.49081	0.49118	0.49155	0.49192	0.49229	0.49266	0.49303	0.4934	0.49377	0.49414	0.49451	0.49488	0.49525	0.49562	0.49599	0.49636	0.49673	0.4971	0.49747	0.49784	0.49821	0.49858	0.49895	0.49933	0.4997	0.50007	0.50044	0.50081	0.50118	0.50155	0.50192	0.50229	0.50266	0.50303	0.5034	0.50377	0.50414	0.50451	0.50488	0.50525	0.50562	0.50599	0.50636	0.50673	0.5071	0.50747	0.50784	0.50821	0.50858	0.50895	0.50932	0.50969	0.51006	0.51043	0.5108	0.51117	0.51154	0.51191	0.51228	0.51265	0.51302	0.51339	0.51376	0.51413	0.5145	0.51487	0.51524	0.51561	0.51598	0.51635	0.51672	0.51709	0.51746	0.51783	0.5182	0.51857	0.51894	0.51931	0.51968	0.52005	0.52042	0.52079	0.52116	0.52153	0.5219	0.52227	0.52264	0.52301	0.52338	0.52375	0.52412	0.52449	0.52486	0.52523	0.5256	0.52597	0.52634	0.52671	0.52708	0.52745	0.52782	0.52819	0.52856	0.52893	0.5293	0.52967	0.53004	0.53041	0.53078	0.53115	0.53152	0.53189	0.53226	0.53263	0.533	0.53337	0.53374	0.53411	0.53448	0.53485	0.53522	0.53559	0.53596	0.53633	0.5367	0.53707	0.53744	0.53781	0.53818	0.53855	0.53892	0.53929	0.53966	0.54003	0.5404	0.54077	0.54114	0.54151	0.54188	0.54225	0.54262	0.54299	0.54336	0.54373	0.5441	0.54447	0.54484	0.54521	0.54558	0.54595	0.54632	0.5467	0.54707	0.54744	0.54781	0.54818	0.54855	0.54892	0.54929	0.54966	0.55003	0.5504	0.55077	0.55114	0.55151	0.55188	0.55225	0.55262	0.55299	0.55336	0.55373	0.5541	0.55447	0.55484	0.55521	0.55558	0.55595	0.55632	0.55669	0.55706	0.55743	0.5578	0.55817	0.55854	0.55891	0.55928	0.55965	0.56002	0.56039	0.56076	0.56113	0.5615	0.56187	0.56224	0.56261	0.56298	0.56335	0.56372	0.56409	0.56446	0.56483	0.5652	0.56557	0.56594	0.56631	0.56668	0.56705	0.56742	0.56779	0.56816	0.56853	0.5689	0.56927	0.56964	0.57001	0.57038	0.57075	0.57112	0.57149	0.57186	0.57223	0.5726	0.57297	0.57334	0.57371	0.57408	0.57445	0.57482	0.57519	0.57556	0.57593	0.5763	0.57667	0.57704	0.57741	0.57778	0.57815	0.57852	0.57889	0.57926	0.57963	0.58	0.58037	0.58074	0.58111	0.58148	0.58185	0.58222	0.58259	0.58296	0.58333	0.5837	0.58407	0.58444	0.58481	0.58518	0.58555	0.58592	0.58629	0.58666	0.58703	0.5874	0.58777	0.58814	0.58851	0.58888	0.58925	0.58962	0.58999	0.59036	0.59073	0.5911	0.59147	0.59184	0.59221	0.59258	0.59295	0.59332	0.59369	0.59407	0.59444	0.59481	0.59518	0.59555	0.59592	0.59629	0.59666	0.59703	0.5974	0.59777	0.59814	0.59851	0.59888	0.59925	0.59962	0.59999	0.60036	0.60073	0.6011	0.60147	0.60184	0.60221	0.60258	0.60295	0.60332	0.60369	0.60406	0.60443	0.6048	0.60517	0.60554	0.60591	0.60628	0.60665	0.60702	0.60739	0.60776	0.60813	0.6085	0.60887	0.60924	0.60961	0.60998	0.61035	0.61072	0.61109	0.61146	0.61183	0.6122	0.61257	0.61294	0.61331	0.61368	0.61405	0.61442	0.61479	0.61516	0.61553	0.6159	0.61627	0.61664	0.61701	0.61738	0.61775	0.61812	0.61849	0.61886	0.61923	0.6196	0.61997	0.62034	0.62071	0.62108	0.62145	0.62182	0.62219	0.62256	0.62293	0.6233	0.62367	0.62404	0.62441	0.62478	0.62515	0.62552	0.62589	0.62626	0.62663	0.627	0.62737	0.62774	0.62811	0.62848	0.62885	0.62922	0.62959	0.62996	0.63033	0.6307	0.63107	0.63144	0.63181	0.63218	0.63255	0.63292	0.63329	0.63366	0.63403	0.6344	0.63477	0.63514	0.63551	0.63588	0.63625	0.63662	0.63699	0.63736	0.63773	0.6381	0.63847	0.63884	0.63921	0.63958	0.63995	0.64032	0.6407	0.64107	0.64144	0.64181	0.64218	0.64255	0.64292	0.64329	0.64366	0.64403	0.6444	0.64477	0.64514	0.64551	0.64588	0.64625	0.64662	0.64699	0.64736	0.64773	0.6481	0.64847	0.64884	0.64921	0.64958	0.64995	0.65032	0.65069	0.65106	0.65143	0.6518	0.65217	0.65254	0.65291	0.65328	0.65365	0.65402	0.65439	0.65476	0.65513	0.6555	0.65587	0.65624	0.65661	0.65698	0.65735	0.65772	0.65809	0.65846	0.65883	0.6592	0.65957	0.65994	0.66031	0.66068	0.66105	0.66142	0.66179	0.66216	0.66253	0.6629	0.66327	0.66364	0.66401	0.66438	0.66475	0.66512	0.66549	0.66586	0.66623	0.6666	0.66697	0.66734	0.66771	0.66808	0.66845	0.66882	0.66919	0.66956	0.66993	0.6703	0.67067	0.67104	0.67141	0.67178	0.67215	0.67252	0.67289	0.67326	0.67363	0.674	0.67437	0.67474	0.67511	0.67548	0.67585	0.67622	0.67659	0.67696	0.67733	0.6777	0.67807	0.67844	0.67881	0.67918	0.67955	0.67992	0.68029	0.68066	0.68103	0.6814	0.68177	0.68214	0.68251	0.68288	0.68325	0.68362	0.68399	0.68436	0.68473	0.6851	0.68547	0.68584	0.68621	0.68658	0.68695	0.68732	0.68769	0.68807	0.68844	0.68881	0.68918	0.68955	0.68992	0.69029	0.69066	0.69103	0.6914	0.69177	0.69214	0.69251	0.69288	0.69325	0.69362	0.69399	0.69436	0.69473	0.6951	0.69547	0.69584	0.69621	0.69658	0.69695	0.69732	0.69769	0.69806	0.69843	0.6988	0.69917	0.69954	0.69991	0.70028	0.70065	0.70102	0.70139	0.70176	0.70213	0.7025	0.70287	0.70324	0.70361	0.70398	0.70435	0.70472	0.70509	0.70546	0.70583	0.7062	0.70657	0.70694	0.70731	0.70768	0.70805	0.70842	0.70879	0.70916	0.70953	0.7099	0.71027	0.71064	0.71101	0.71138	0.71175	0.71212	0.71249	0.71286	0.71323	0.7136	0.71397	0.71434	0.71471	0.71508	0.71545	0.71582	0.71619	0.71656	0.71693	0.7173	0.71767	0.71804	0.71841	0.71878	0.71915	0.71952	0.71989	0.72026	0.72063	0.721	0.72137	0.72174	0.72211	0.72248	0.72285	0.72322	0.72359	0.72396	0.72433	0.7247	0.72507	0.72544	0.72581	0.72618	0.72655	0.72692	0.72729	0.72766	0.72803	0.7284	0.72877	0.72914	0.72951	0.72988	0.73025	0.73062	0.73099	0.73136	0.73173	0.7321	0.73247	0.73284	0.73321	0.73358	0.73395	0.73432	0.73469	0.73506	0.73544	0.73581	0.73618	0.73655	0.73692	0.73729	0.73766	0.73803	0.7384	0.73877	0.73914	0.73951	0.73988	0.74025	0.74062	0.74099	0.74136	0.74173	0.7421	0.74247	0.74284	0.74321	0.74358	0.74395	0.74432	0.74469	0.74506	0.74543	0.7458	0.74617	0.74654	0.74691	0.74728	0.74765	0.74802	0.74839	0.74876	0.74913	0.7495	0.74987	0.75024	0.75061	0.75098	0.75135	0.75172	0.75209	0.75246	0.75283	0.7532	0.75357	0.75394	0.75431	0.75468	0.75505	0.75542	0.75579	0.75616	0.75653	0.7569	0.75727	0.75764	0.75801	0.75838	0.75875	0.75912	0.75949	0.75986	0.76023	0.7606	0.76097	0.76134	0.76171	0.76208	0.76245	0.76282	0.76319	0.76356	0.76393	0.7643	0.76467	0.76504	0.76541	0.76578	0.76615	0.76652	0.76689	0.76726	0.76763	0.768	0.76837	0.76874	0.76911	0.76948	0.76985	0.77022	0.77059	0.77096	0.77133	0.7717	0.77207	0.77244	0.77281	0.77318	0.77355	0.77392	0.77429	0.77466	0.77503	0.7754	0.77577	0.77614	0.77651	0.77688	0.77725	0.77762	0.77799	0.77836	0.77873	0.7791	0.77947	0.77984	0.78021	0.78058	0.78095	0.78132	0.78169	0.78207	0.78243	0.78281	0.78318	0.78355	0.78392	0.78429	0.78466	0.78503	0.7854	0.78577	0.78614	0.78651	0.78688	0.78725	0.78762	0.78799	0.78836	0.78873	0.7891	0.78947	0.78984	0.79021	0.79058	0.79095	0.79132	0.79169	0.79206	0.79243	0.7928	0.79317	0.79354]*1e-6;
n=[3.4144	3.4154	3.4165	3.4177	3.419	3.4203	3.4218	3.4234	3.4251	3.4269	3.4289	3.431	3.4332	3.4356	3.4382	3.4409	3.4438	3.4468	3.45	3.4534	3.457	3.4608	3.4648	3.469	3.4734	3.478	3.4827	3.4878	3.493	3.4984	3.504	3.5098	3.5158	3.5221	3.5285	3.5351	3.5418	3.5487	3.5558	3.563	3.5704	3.5779	3.5854	3.5931	3.6008	3.6086	3.6165	3.6243	3.6322	3.64	3.6478	3.6555	3.6632	3.6708	3.6782	3.6855	3.6926	3.6996	3.7064	3.7129	3.7192	3.7252	3.731	3.7365	3.7417	3.7466	3.7512	3.7554	3.7593	3.7628	3.766	3.7689	3.7713	3.7735	3.7753	3.7767	3.7778	3.7785	3.7789	3.7791	3.7789	3.7784	3.7776	3.7766	3.7753	3.7738	3.772	3.7701	3.768	3.7658	3.7634	3.7609	3.7583	3.7557	3.753	3.7503	3.7475	3.7448	3.7421	3.7394	3.7368	3.7343	3.7318	3.7295	3.7273	3.7252	3.7233	3.7215	3.7198	3.7184	3.7171	3.716	3.7151	3.7144	3.7138	3.7135	3.7133	3.7134	3.7136	3.714	3.7146	3.7154	3.7164	3.7175	3.7188	3.7202	3.7219	3.7236	3.7255	3.7275	3.7297	3.732	3.7343	3.7368	3.7394	3.742	3.7448	3.7476	3.7504	3.7533	3.7563	3.7593	3.7623	3.7654	3.7684	3.7715	3.7746	3.7776	3.7807	3.7837	3.7868	3.7898	3.7927	3.7956	3.7985	3.8013	3.8041	3.8068	3.8095	3.8121	3.8146	3.8171	3.8195	3.8218	3.824	3.8261	3.8282	3.8302	3.832	3.8338	3.8355	3.8371	3.8386	3.84	3.8413	3.8425	3.8435	3.8445	3.8454	3.8461	3.8468	3.8473	3.8477	3.848	3.8482	3.8483	3.8483	3.8481	3.8479	3.8475	3.847	3.8464	3.8457	3.8449	3.8439	3.8429	3.8418	3.8406	3.8394	3.8382	3.837	3.8359	3.835	3.8345	3.8344	3.8349	3.8364	3.839	3.843	3.8488	3.8566	3.8668	3.8795	3.895	3.9131	3.9338	3.9567	3.9814	4.0069	4.0325	4.0571	4.0796	4.0988	4.1137	4.1234	4.1272	4.1247	4.1159	4.1009	4.0804	4.0551	4.0261	3.9944	3.9612	3.9276	3.8945	3.8629	3.8333	3.8062	3.7818	3.7601	3.7412	3.7248	3.7107	3.6985	3.688	3.6788	3.6708	3.6636	3.657	3.6509	3.6452	3.6397	3.6344	3.6292	3.6241	3.619	3.614	3.609	3.604	3.5991	3.5941	3.5892	3.5843	3.5794	3.5745	3.5697	3.5649	3.56	3.5553	3.5505	3.5458	3.5411	3.5364	3.5317	3.5271	3.5226	3.518	3.5135	3.509	3.5046	3.5002	3.4959	3.4916	3.4873	3.4831	3.4789	3.4748	3.4707	3.4666	3.4626	3.4587	3.4548	3.4509	3.4471	3.4434	3.4397	3.436	3.4324	3.4288	3.4253	3.4218	3.4184	3.4151	3.4117	3.4085	3.4053	3.4021	3.399	3.3959	3.3929	3.3899	3.387	3.3841	3.3813	3.3785	3.3758	3.3731	3.3705	3.3679	3.3653	3.3628	3.3604	3.358	3.3556	3.3533	3.351	3.3488	3.3466	3.3444	3.3423	3.3402	3.3382	3.3362	3.3342	3.3323	3.3304	3.3286	3.3267	3.325	3.3232	3.3215	3.3198	3.3182	3.3166	3.315	3.3134	3.3119	3.3104	3.3089	3.3075	3.3061	3.3047	3.3033	3.302	3.3007	3.2994	3.2981	3.2969	3.2957	3.2945	3.2933	3.2921	3.291	3.2899	3.2888	3.2877	3.2866	3.2856	3.2846	3.2835	3.2825	3.2816	3.2806	3.2797	3.2787	3.2778	3.2769	3.276	3.2751	3.2742	3.2733	3.2725	3.2717	3.2708	3.27	3.2692	3.2684	3.2676	3.2668	3.266	3.2652	3.2645	3.2637	3.263	3.2622	3.2615	3.2608	3.26	3.2593	3.2586	3.2579	3.2572	3.2565	3.2558	3.2551	3.2544	3.2537	3.2531	3.2524	3.2517	3.251	3.2504	3.2497	3.249	3.2484	3.2477	3.2471	3.2464	3.2458	3.2451	3.2445	3.2438	3.2432	3.2425	3.2419	3.2412	3.2406	3.24	3.2393	3.2387	3.2381	3.2374	3.2368	3.2361	3.2355	3.2349	3.2342	3.2336	3.233	3.2323	3.2317	3.2311	3.2304	3.2298	3.2292	3.2285	3.2279	3.2273	3.2266	3.226	3.2254	3.2247	3.2241	3.2235	3.2228	3.2222	3.2215	3.2209	3.2203	3.2196	3.219	3.2183	3.2177	3.2171	3.2164	3.2158	3.2151	3.2145	3.2138	3.2132	3.2125	3.2119	3.2112	3.2106	3.2099	3.2093	3.2086	3.208	3.2073	3.2067	3.206	3.2054	3.2047	3.204	3.2034	3.2027	3.2021	3.2014	3.2007	3.2001	3.1994	3.1987	3.1981	3.1974	3.1967	3.1961	3.1954	3.1947	3.1941	3.1934	3.1927	3.192	3.1914	3.1907	3.19	3.1893	3.1887	3.188	3.1873	3.1866	3.1859	3.1852	3.1846	3.1839	3.1832	3.1825	3.1818	3.1811	3.1804	3.1797	3.179	3.1784	3.1777	3.177	3.1763	3.1756	3.1749	3.1742	3.1735	3.1728	3.1721	3.1714	3.1707	3.17	3.1693	3.1686	3.1678	3.1671	3.1664	3.1657	3.165	3.1643	3.1636	3.1629	3.1622	3.1614	3.1607	3.16	3.1593	3.1586	3.1578	3.1571	3.1564	3.1557	3.155	3.1542	3.1535	3.1528	3.152	3.1513	3.1506	3.1499	3.1491	3.1484	3.1477	3.1469	3.1462	3.1455	3.1447	3.144	3.1432	3.1425	3.1418	3.141	3.1403	3.1395	3.1388	3.1381	3.1373	3.1366	3.1358	3.1351	3.1343	3.1336	3.1328	3.1321	3.1313	3.1306	3.1298	3.129	3.1283	3.1275	3.1268	3.126	3.1253	3.1245	3.1237	3.123	3.1222	3.1214	3.1207	3.1199	3.1191	3.1184	3.1176	3.1168	3.1161	3.1153	3.1145	3.1137	3.113	3.1122	3.1114	3.1106	3.1099	3.1091	3.1083	3.1075	3.1067	3.106	3.1052	3.1044	3.1036	3.1028	3.102	3.1013	3.1005	3.0997	3.0989	3.0981	3.0973	3.0965	3.0957	3.0949	3.0941	3.0933	3.0925	3.0917	3.0909	3.0901	3.0893	3.0885	3.0877	3.0869	3.0861	3.0853	3.0845	3.0837	3.0829	3.0821	3.0813	3.0805	3.0797	3.0789	3.078	3.0772	3.0764	3.0756	3.0748	3.074	3.0732	3.0723	3.0715	3.0707	3.0699	3.0691	3.0682	3.0674	3.0666	3.0658	3.0649	3.0641	3.0633	3.0625	3.0616	3.0608	3.06	3.0591	3.0583	3.0575	3.0566	3.0558	3.055	3.0541	3.0533	3.0525	3.0516	3.0508	3.0499	3.0491	3.0483	3.0474	3.0466	3.0457	3.0449	3.044	3.0432	3.0423	3.0415	3.0406	3.0398	3.0389	3.0381	3.0372	3.0364	3.0355	3.0347	3.0338	3.033	3.0321	3.0312	3.0304	3.0295	3.0287	3.0278	3.0269	3.0261	3.0252	3.0244	3.0235	3.0226	3.0218	3.0209	3.02	3.0192	3.0183	3.0174	3.0165	3.0157	3.0148	3.0139	3.013	3.0122	3.0113	3.0104	3.0095	3.0087	3.0078	3.0069	3.006	3.0051	3.0043	3.0034	3.0025	3.0016	3.0007	2.9998	2.9989	2.9981	2.9972	2.9963	2.9954	2.9945	2.9936	2.9927	2.9918	2.9909	2.99	2.9891	2.9882	2.9873	2.9864	2.9855	2.9846	2.9837	2.9828	2.9819	2.981	2.9801	2.9792	2.9783	2.9774	2.9765	2.9756	2.9747	2.9738	2.9729	2.972	2.9711	2.9702	2.9692	2.9683	2.9674	2.9665	2.9656	2.9647	2.9638	2.9628	2.9619	2.961	2.9601	2.9592	2.9582	2.9573	2.9564	2.9555	2.9545	2.9536	2.9527	2.9518	2.9508	2.9499	2.949	2.9481	2.9471	2.9462	2.9453	2.9443	2.9434	2.9425	2.9415	2.9406	2.9397	2.9387	2.9378	2.9369	2.9359	2.935	2.934	2.9331	2.9322	2.9312	2.9303	2.9293	2.9284	2.9274	2.9265	2.9256	2.9246	2.9237	2.9227	2.9218	2.9208	2.9199	2.9189	2.918	2.917	2.9161	2.9151	2.9142	2.9132	2.9122	2.9113	2.9103	2.9094	2.9084	2.9075	2.9065	2.9055	2.9046	2.9036	2.9026	2.9017	2.9007	2.8998	2.8988	2.8978	2.8969	2.8959	2.8949	2.894	2.893	2.892	2.891	2.8901	2.8891	2.8881	2.8872	2.8862	2.8852	2.8842	2.8833	2.8823	2.8813	2.8803	2.8793	2.8784	2.8774	2.8764	2.8754	2.8744	2.8735	2.8725	2.8715	2.8705	2.8695	2.8685	2.8676	2.8666	2.8656	2.8646	2.8636	2.8626	2.8616	2.8606	2.8596	2.8586	2.8577	2.8567	2.8557	2.8547	2.8537	2.8527	2.8517	2.8507	2.8497	2.8487	2.8477	2.8467	2.8457	2.8447	2.8437	2.8427	2.8417	2.8407	2.8397	2.8387	2.8377	2.8367	2.8357	2.8346	2.8336	2.8326	2.8316	2.8306	2.8296	2.8286	2.8276	2.8266	2.8256	2.8245	2.8235	2.8225	2.8215	2.8205	2.8195	2.8185	2.8174	2.8164	2.8154	2.8144	2.8134	2.8123	2.8113	2.8103	2.8093	2.8082	2.8072	2.8062	2.8052	2.8042	2.8031	2.8021	2.8011	2.8	2.799	2.798	2.797	2.7959	2.7949	2.7939	2.7928	2.7918	2.7908	2.7897	2.7887	2.7877	2.7866	2.7856	2.7846	2.7835	2.7825	2.7814	2.7804	2.7794	2.7783	2.7773	2.7762	2.7752	2.7742	2.7731	2.7721	2.771	2.77	2.7689	2.7679	2.7669	2.7658	2.7648	2.7637	2.7627	2.7616	2.7606	2.7595	2.7585	2.7574	2.7564	2.7553	2.7543	2.7532];
k=[1.6123	1.6019	1.5915	1.581	1.5705	1.56	1.5494	1.5388	1.5282	1.5175	1.5069	1.4962	1.4855	1.4748	1.464	1.4533	1.4425	1.4317	1.4209	1.4101	1.3993	1.3885	1.3777	1.3669	1.356	1.3452	1.3344	1.3235	1.3127	1.3019	1.2911	1.2803	1.2695	1.2587	1.248	1.2372	1.2265	1.2157	1.205	1.1943	1.1836	1.173	1.1623	1.1517	1.1411	1.1306	1.12	1.1095	1.099	1.0886	1.0782	1.0678	1.0574	1.0471	1.0368	1.0265	1.0163	1.0061	0.996	0.9859	0.9758	0.9658	0.9558	0.9458	0.9359	0.9261	0.9163	0.9065	0.8968	0.8871	0.8775	0.8679	0.8584	0.8489	0.8395	0.8301	0.8208	0.8115	0.8023	0.7932	0.7841	0.775	0.766	0.7571	0.7482	0.7393	0.7306	0.7219	0.7132	0.7046	0.6961	0.6876	0.6792	0.6708	0.6625	0.6543	0.6461	0.638	0.6299	0.622	0.614	0.6062	0.5984	0.5906	0.583	0.5754	0.5678	0.5604	0.553	0.5456	0.5383	0.5312	0.524	0.517	0.51	0.5031	0.4963	0.4895	0.4828	0.4763	0.4698	0.4634	0.457	0.4508	0.4447	0.4387	0.4328	0.427	0.4213	0.4158	0.4104	0.4051	0.4	0.3951	0.3903	0.3857	0.3813	0.3772	0.3732	0.3694	0.366	0.3627	0.3598	0.3571	0.3547	0.3527	0.351	0.3496	0.3487	0.3481	0.3479	0.3481	0.3488	0.3499	0.3515	0.3535	0.356	0.359	0.3625	0.3665	0.371	0.3759	0.3814	0.3873	0.3937	0.4006	0.4078	0.4155	0.4236	0.4321	0.4408	0.4499	0.4592	0.4687	0.4783	0.4881	0.4979	0.5078	0.5175	0.5272	0.5367	0.5459	0.5548	0.5634	0.5716	0.5792	0.5864	0.5929	0.5988	0.604	0.6085	0.6121	0.6149	0.6169	0.618	0.6182	0.6174	0.6157	0.6131	0.6095	0.605	0.5996	0.5932	0.586	0.5779	0.569	0.5594	0.5489	0.5378	0.5261	0.5137	0.5009	0.4875	0.4738	0.4597	0.4453	0.4306	0.4158	0.401	0.3861	0.3717	0.359	0.3515	0.3566	0.3854	0.4473	0.5383	0.6316	0.6853	0.6678	0.5815	0.4606	0.3465	0.2631	0.2124	0.1844	0.1681	0.1568	0.1474	0.1388	0.1307	0.1231	0.116	0.1092	0.1029	0.0969	0.0913	0.0861	0.0812	0.0767	0.0725	0.0685	0.0648	0.0614	0.0583	0.0553	0.0526	0.0501	0.0477	0.0455	0.0435	0.0416	0.0399	0.0383	0.0367	0.0353	0.034	0.0328	0.0316	0.0305	0.0295	0.0286	0.0277	0.0268	0.026	0.0252	0.0245	0.0238	0.0231	0.0225	0.0219	0.0213	0.0207	0.0202	0.0196	0.0191	0.0186	0.0181	0.0177	0.0172	0.0168	0.0164	0.0159	0.0155	0.0151	0.0148	0.0144	0.014	0.0137	0.0133	0.013	0.0126	0.0123	0.012	0.0117	0.0114	0.0111	0.0108	0.0105	0.0103	0.01	0.0097	0.0095	0.0092	0.009	0.0087	0.0085	0.0083	0.0081	0.0079	0.0076	0.0074	0.0072	0.0071	0.0069	0.0067	0.0065	0.0063	0.0062	0.006	0.0058	0.0057	0.0055	0.0054	0.0052	0.0051	0.0049	0.0048	0.0047	0.0045	0.0044	0.0043	0.0042	0.004	0.0039	0.0038	0.0037	0.0036	0.0035	0.0034	0.0033	0.0032	0.0031	0.003	0.0029	0.0029	0.0028	0.0027	0.0026	0.0025	0.0025	0.0024	0.0023	0.0023	0.0022	0.0021	0.0021	0.002	0.0019	0.0019	0.0018	0.0018	0.0017	0.0017	0.0016	0.0016	0.0015	0.0015	0.0014	0.0014	0.0013	0.0013	0.0013	0.0012	0.0012	0.0012	0.0011	0.0011	0.0011	1E-3	1E-3	1E-3	9E-4	9E-4	9E-4	8E-4	8E-4	8E-4	8E-4	7E-4	7E-4	7E-4	7E-4	7E-4	6E-4	6E-4	6E-4	6E-4	6E-4	5E-4	5E-4	5E-4	5E-4	5E-4	5E-4	4E-4	4E-4	4E-4	4E-4	4E-4	4E-4	4E-4	3E-4	3E-4	3E-4	3E-4	3E-4	3E-4	3E-4	3E-4	3E-4	3E-4	2E-4	2E-4	2E-4	2E-4	2E-4	2E-4	2E-4	2E-4	2E-4	2E-4	2E-4	2E-4	2E-4	2E-4	2E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	1E-4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0];

aux=interp1(longdaexp,n,lambda);
aux2=interp1(longdaexp,k,lambda);

N=aux+aux2*1i;
end


function N=EVA(lambda)
%indice de refraccion del eva
%Conventional ethylene-vinyl acetate (EVA) with high UV absorption. Room temperature.
%M. R. Vogt, H. Schulte-Huxel, D. Hinken, H. Holst, M. Winter, S. Blankemeyer, R. Witteck, M. Köntges, K. Bothe, R. Brendel. Optical constants of UV transparent EVA and the impact on the PV module output power under realistic illumination, Energy Procedia 92, 523-530 (2016)

longdaexp=[0.25	0.26	0.27	0.28	0.29	0.3	0.31	0.32	0.33	0.34	0.35	0.36	0.37	0.38	0.39	0.4	0.41	0.42	0.43	0.44	0.45	0.46	0.47	0.48	0.49	0.5	0.51	0.52	0.53	0.54	0.55	0.56	0.57	0.58	0.59	0.6	0.61	0.62	0.63	0.64	0.65	0.66	0.67	0.68	0.69	0.7	0.71	0.72	0.73	0.74	0.75	0.76	0.77	0.78	0.79	0.8	0.81	0.82	0.83	0.84	0.85	0.86	0.87	0.88	0.89	0.9	0.91	0.92	0.93	0.94	0.95	0.96	0.97	0.98	0.99	1	1.01	1.02	1.03	1.04	1.05	1.06	1.07	1.08	1.09	1.1	1.11	1.12	1.13	1.14	1.15	1.16	1.17	1.18	1.19	1.2	1.21	1.22	1.23	1.24	1.25	1.26	1.27	1.28	1.29	1.3	1.31	1.32	1.33	1.34	1.35	1.36	1.37	1.38	1.39	1.4	1.41	1.42	1.43	1.44	1.45	1.46	1.47	1.48	1.49	1.5	1.51	1.52	1.53	1.54	1.55	1.56	1.57	1.58	1.59	1.6	1.61	1.62	1.63	1.64	1.65	1.66	1.67	1.68	1.69	1.7	1.71	1.72	1.73	1.74	1.75	1.76	1.77	1.78	1.79	1.8	1.81	1.82	1.83	1.84	1.85	1.86	1.87	1.88	1.89	1.9	1.91	1.92	1.93	1.94	1.95	1.96	1.97	1.98	1.99	2	2.01	2.02	2.03	2.04	2.05	2.06	2.07	2.08	2.09	2.1	2.11	2.12	2.13	2.14	2.15	2.16	2.17	2.18	2.19	2.2	2.21	2.22	2.23	2.24	2.25	2.26	2.27	2.28	2.29	2.3	2.31	2.32	2.33	2.34	2.35	2.36	2.37	2.38	2.39	2.4	2.41	2.42	2.43	2.44	2.45	2.46	2.47	2.48	2.49	2.5]*1e-6;
n=[1.562	1.555	1.549	1.543	1.538	1.533	1.529	1.525	1.522	1.519	1.517	1.514	1.512	1.51	1.508	1.507	1.505	1.504	1.503	1.502	1.501	1.5	1.499	1.498	1.497	1.496	1.496	1.495	1.494	1.494	1.493	1.493	1.492	1.492	1.491	1.491	1.49	1.49	1.49	1.489	1.489	1.489	1.488	1.488	1.488	1.488	1.487	1.487	1.487	1.487	1.487	1.486	1.486	1.486	1.486	1.486	1.486	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.482	1.482	1.482	1.482	1.482	1.482	1.482	1.482	1.482	1.482	1.482	1.482	1.482	1.482	1.482	1.482	1.482	1.482	1.482	1.482	1.482	1.482	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.481	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48	1.48];
k=[7.71E-4	0.00189	5.92E-4	0.00142	0.00216	0.00314	0.0034	0.00303	0.00105	0.00456	0.00187	1.1E-4	4.45E-5	1.65E-5	6.61E-6	2.62E-6	1.65E-6	1.28E-6	1.11E-6	1.02E-6	9.66E-7	9.2E-7	8.8E-7	8.43E-7	8.11E-7	7.79E-7	7.49E-7	7.24E-7	6.98E-7	6.74E-7	6.52E-7	6.29E-7	6.1E-7	5.89E-7	5.71E-7	5.52E-7	5.37E-7	5.19E-7	5.03E-7	4.91E-7	4.79E-7	4.6E-7	4.42E-7	4.31E-7	4.18E-7	4.08E-7	3.98E-7	3.9E-7	3.85E-7	3.88E-7	3.97E-7	3.99E-7	3.78E-7	3.47E-7	3.27E-7	3.18E-7	3.18E-7	3.19E-7	3.21E-7	3.18E-7	3.08E-7	3E-7	3.05E-7	3.3E-7	3.91E-7	4.84E-7	6.49E-7	9.1E-7	1.1E-6	9.66E-7	6.1E-7	3.51E-7	2.98E-7	2.84E-7	3.06E-7	3.71E-7	4.59E-7	5.51E-7	6.48E-7	7.03E-7	6.73E-7	5.79E-7	4.9E-7	4.19E-7	3.77E-7	3.65E-7	3.93E-7	4.77E-7	6.69E-7	1.12E-6	1.97E-6	3.22E-6	4.66E-6	6.33E-6	8.87E-6	1.24E-5	1.51E-5	1.27E-5	8.18E-6	4.87E-6	3.18E-6	2.27E-6	1.78E-6	1.51E-6	1.33E-6	1.2E-6	1.08E-6	9.99E-7	1.04E-6	1.3E-6	1.75E-6	2.63E-6	4.36E-6	7.13E-6	9.89E-6	1.1E-5	1.14E-5	1.16E-5	1.1E-5	9.69E-6	8.24E-6	6.84E-6	5.56E-6	4.49E-6	3.74E-6	3.27E-6	3.06E-6	3.06E-6	3.15E-6	3.1E-6	2.74E-6	2.39E-6	2.17E-6	2.08E-6	2.1E-6	2.23E-6	2.47E-6	2.9E-6	3.65E-6	4.77E-6	6.55E-6	9.31E-6	1.33E-5	1.82E-5	2.39E-5	3.34E-5	5.16E-5	8.7E-5	9.88E-5	8.06E-5	7.82E-5	7.96E-5	6.85E-5	5.82E-5	5.42E-5	5.34E-5	5.29E-5	5.07E-5	4.74E-5	4.33E-5	3.92E-5	3.58E-5	3.33E-5	3.22E-5	3.27E-5	3.32E-5	3.33E-5	3.33E-5	3.38E-5	3.45E-5	3.37E-5	3.27E-5	3.23E-5	3.23E-5	3.3E-5	3.44E-5	3.49E-5	3.37E-5	3.27E-5	3.28E-5	3.28E-5	3.28E-5	3.22E-5	3.22E-5	3.37E-5	3.53E-5	3.71E-5	3.97E-5	4.48E-5	4.82E-5	4.32E-5	3.92E-5	3.8E-5	3.86E-5	3.95E-5	4.21E-5	4.76E-5	5.55E-5	6.91E-5	8.89E-5	1.57E-4	2.11E-4	2.7E-4	3.59E-4	5.18E-4	7.34E-4	8.12E-4	6.54E-4	5.95E-4	6.69E-4	7.58E-4	7.35E-4	7.06E-4	7.47E-4	7.56E-4	7.5E-4	7.14E-4	6.84E-4	6.63E-4	6.21E-4	5.78E-4	5.33E-4	5.03E-4	4.65E-4	4.45E-4	4.19E-4];

aux=interp1(longdaexp,n,lambda);
aux2=interp1(longdaexp,k,lambda);

N=aux+aux2*1i;

end

function N=EVA_UV(lambda)
%indice de refraccion del eva
%Ethylene-vinyl acetate (EVA) with enhanced UV transmission. Room temperature.
%M. R. Vogt, H. Schulte-Huxel, D. Hinken, H. Holst, M. Winter, S. Blankemeyer, R. Witteck, M. Köntges, K. Bothe, R. Brendel. Optical constants of UV transparent EVA and the impact on the PV module output power under realistic illumination, Energy Procedia 92, 523-530 (2016)

longdaexp=[0.25	0.26	0.27	0.28	0.29	0.3	0.31	0.32	0.33	0.34	0.35	0.36	0.37	0.38	0.39	0.4	0.41	0.42	0.43	0.44	0.45	0.46	0.47	0.48	0.49	0.5	0.51	0.52	0.53	0.54	0.55	0.56	0.57	0.58	0.59	0.6	0.61	0.62	0.63	0.64	0.65	0.66	0.67	0.68	0.69	0.7	0.71	0.72	0.73	0.74	0.75	0.76	0.77	0.78	0.79	0.8	0.81	0.82	0.83	0.84	0.85	0.86	0.87	0.88	0.89	0.9	0.91	0.92	0.93	0.94	0.95	0.96	0.97	0.98	0.99	1	1.01	1.02	1.03	1.04	1.05	1.06	1.07	1.08	1.09	1.1	1.11	1.12	1.13	1.14	1.15	1.16	1.17	1.18	1.19	1.2	1.21	1.22	1.23	1.24	1.25	1.26	1.27	1.28	1.29	1.3	1.31	1.32	1.33	1.34	1.35	1.36	1.37	1.38	1.39	1.4	1.41	1.42	1.43	1.44	1.45	1.46	1.47	1.48	1.49	1.5	1.51	1.52	1.53	1.54	1.55	1.56	1.57	1.58	1.59	1.6	1.61	1.62	1.63	1.64	1.65	1.66	1.67	1.68	1.69	1.7	1.71	1.72	1.73	1.74	1.75	1.76	1.77	1.78	1.79	1.8	1.81	1.82	1.83	1.84	1.85	1.86	1.87	1.88	1.89	1.9	1.91	1.92	1.93	1.94	1.95	1.96	1.97	1.98	1.99	2	2.01	2.02	2.03	2.04	2.05	2.06	2.07	2.08	2.09	2.1	2.11	2.12	2.13	2.14	2.15	2.16	2.17	2.18	2.19	2.2	2.21	2.22	2.23	2.24	2.25	2.26	2.27	2.28	2.29	2.3	2.31	2.32	2.33	2.34	2.35	2.36	2.37	2.38	2.39	2.4	2.41	2.42	2.43	2.44	2.45	2.46	2.47	2.48	2.49	2.5]*1e-6;
n=[1.557	1.551	1.545	1.54	1.536	1.532	1.528	1.525	1.522	1.519	1.517	1.515	1.513	1.511	1.51	1.508	1.507	1.506	1.505	1.504	1.503	1.502	1.501	1.5	1.499	1.499	1.498	1.497	1.497	1.496	1.496	1.495	1.495	1.494	1.494	1.493	1.493	1.493	1.492	1.492	1.492	1.491	1.491	1.491	1.491	1.49	1.49	1.49	1.49	1.49	1.489	1.489	1.489	1.489	1.489	1.488	1.488	1.488	1.488	1.488	1.488	1.488	1.487	1.487	1.487	1.487	1.487	1.487	1.487	1.487	1.487	1.487	1.486	1.486	1.486	1.486	1.486	1.486	1.486	1.486	1.486	1.486	1.486	1.486	1.486	1.486	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.485	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.484	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483	1.483];
k=[9.86E-6	9.7E-6	8.2E-6	6.64E-6	5.31E-6	3.89E-6	3.11E-6	2.64E-6	2.28E-6	1.95E-6	1.56E-6	1.27E-6	1.07E-6	9.48E-7	8.57E-7	7.84E-7	7.37E-7	6.96E-7	6.62E-7	6.33E-7	6.04E-7	5.78E-7	5.55E-7	5.32E-7	5.1E-7	4.9E-7	4.7E-7	4.52E-7	4.33E-7	4.17E-7	4.02E-7	3.87E-7	3.74E-7	3.58E-7	3.47E-7	3.34E-7	3.23E-7	3.12E-7	3.03E-7	2.96E-7	2.88E-7	2.74E-7	2.64E-7	2.56E-7	2.49E-7	2.43E-7	2.36E-7	2.32E-7	2.32E-7	2.38E-7	2.53E-7	2.59E-7	2.4E-7	2.1E-7	1.94E-7	1.89E-7	1.91E-7	1.97E-7	2.01E-7	1.98E-7	1.91E-7	1.88E-7	1.91E-7	2.19E-7	2.89E-7	3.86E-7	5.64E-7	8.43E-7	1.05E-6	9.3E-7	5.35E-7	2.72E-7	2.16E-7	2.05E-7	2.35E-7	3E-7	3.98E-7	4.99E-7	6.01E-7	6.68E-7	6.33E-7	5.33E-7	4.4E-7	3.69E-7	3.26E-7	3.11E-7	3.41E-7	4.28E-7	6.35E-7	1.08E-6	1.99E-6	3.31E-6	4.83E-6	6.59E-6	9.23E-6	1.3E-5	1.59E-5	1.34E-5	8.47E-6	5.03E-6	3.26E-6	2.32E-6	1.82E-6	1.53E-6	1.35E-6	1.21E-6	1.08E-6	1E-6	1.04E-6	1.26E-6	1.79E-6	2.66E-6	4.44E-6	7.21E-6	1.04E-5	1.16E-5	1.2E-5	1.23E-5	1.15E-5	1.02E-5	8.63E-6	7.16E-6	5.8E-6	4.66E-6	3.89E-6	3.4E-6	3.19E-6	3.18E-6	3.29E-6	3.24E-6	2.86E-6	2.47E-6	2.25E-6	2.15E-6	2.18E-6	2.3E-6	2.56E-6	3.01E-6	3.78E-6	4.94E-6	6.79E-6	9.65E-6	1.39E-5	1.89E-5	2.5E-5	3.46E-5	5.29E-5	8.61E-5	9.8E-5	8.69E-5	8.12E-5	8.2E-5	7.23E-5	6.15E-5	5.72E-5	5.64E-5	5.58E-5	5.32E-5	5E-5	4.56E-5	4.13E-5	3.78E-5	3.5E-5	3.37E-5	3.41E-5	3.46E-5	3.46E-5	3.47E-5	3.55E-5	3.62E-5	3.55E-5	3.45E-5	3.39E-5	3.39E-5	3.47E-5	3.61E-5	3.66E-5	3.57E-5	3.46E-5	3.45E-5	3.46E-5	3.43E-5	3.38E-5	3.38E-5	3.53E-5	3.74E-5	3.9E-5	4.15E-5	4.72E-5	5.1E-5	4.54E-5	4.1E-5	4.01E-5	4.02E-5	4.13E-5	4.4E-5	4.99E-5	5.81E-5	7.19E-5	9.01E-5	1.57E-4	2.11E-4	2.7E-4	3.59E-4	5.18E-4	7.34E-4	8.12E-4	6.54E-4	5.95E-4	6.69E-4	7.58E-4	7.35E-4	7.06E-4	7.47E-4	7.56E-4	7.5E-4	7.14E-4	6.84E-4	6.63E-4	6.21E-4	5.78E-4	5.33E-4	5.03E-4	4.65E-4	4.45E-4	4.19E-4];

aux=interp1(longdaexp,n,lambda);
aux2=interp1(longdaexp,k,lambda);

N=aux+aux2*1i;

end




%funciones para el calculo de medios efectivos a partir de la combinacion
%de direfentes sustancias

function neff=looyenga(n1,n2,x2)
%esta función calcula el indice de refracción efectivo a partir de los
%índices de refracción de los componentes individuales y la fracción del
%componente 2 utilizando el modelo de looyenga (Pag 117 W. Theis).

e1=n1.^2;
e2=n2.^2;

eeff=(e1.^(1/3).*(1-x2)+e2.^(1/3)*x2).^3;

neff=eeff.^0.5;

end

function neff=looyenga3m(n1,n2,n3,x2,x3)
%esta funciï¿½n calcula el indice de refracciï¿½n efectivo a partir de los
%ï¿½ndices de refracciï¿½n de los componentes individuales y la fracciï¿½n del
%componente 2 y 3 utilizando el modelo de looyenga (Pag 117 W. Theis).

e1=n1.^2;
e2=n2.^2;
e3=n3.^2;

eeff=(e1.^(1/3).*(1-x2-x3)+e2.^(1/3)*x2+e3.^(1/3)*x3).^3;

neff=eeff.^0.5;

end

function neff=looyengacil(n1,n2,x2)
%esta función calcula el indice de refracción efectivo a partir de los
%índices de refracción de los componentes individuales y la fracción del
%componente 2 utilizando el modelo de looyenga (Pag 117 W. Theis) en el
%caso de tener poros cilindricos.

neff=n1*(1-x2)+n2*x2;

end

function neff=mezclahomogenea(f1,n1,n2)
%funcion que devuelve el indice de refraccion de una mezcla de líquidos
%teniendo en cuenta la fraccion en volumen de cada componente y sus indices
%de refracción. Esto corresponde al modelo de Lorentz-Lorenz
%f1 es la fraccion en volumen del componente 1

aux=f1*(n1.^2 - 1)./(n1.^2 + 2) + (1-f1)*(n2.^2 - 1)./(n2.^2 + 2);

neff=sqrt(-1-2*aux)./sqrt(aux-1);
end

function neff=maxwell(n1,n2,x2)
%esta función calcula el indice de refracción efectivo a partir de los
%índices de refracción de los componentes individuales y la fracción del
%componente 2 utilizando el modelo de Maxwell Garnett (Pag 117 W. Theis).

e1=n1.^2;
e2=n2.^2;

eeff=(-3*e1.*e2+2*e1.*e2*x2-2*power(e2,2)*x2)./(-3*e2-e1*x2+e2*x2);
neff=eeff.^0.5;
end

function neff=Bruggeman(n1,n2,x2)
%esta función calcula el indice de refracción efectivo a partir de los
%índices de refracción de los componentes individuales y la fracción del
%componente 2 utilizando el modelo de Bruggeman(Pag 117 W. Theis).

e1=n1.^2;
e2=n2.^2;

eeff=(2*e1-e2-3*e1*x2+3*e2*x2+sqrt(8*e1.*e2+power(2*e1-e2-3*e1*x2+3*e2*x2,2)))/4;

neff=eeff.^0.5;
end
