% =========================================================================
% == Inicio limpio de terminal ============================================
% =========================================================================
clc
path(pathdef)
clear all
close all

addpath(genpath('.'));
% =========================================================================
% == Definicion de la Estructura ==========================================
% =========================================================================
%
% A seguir se define la estructura, se deben utilizar UNIDADES COMPATIBLES.
% NO DEBEN modificarse los nombres de los campos de la variable ES.

% Nombre del problema
ES.Problema = 'Puente';

ES.gamma=1e-3; % Rigidez de los vacios respecto al material estructural.
               % Coincide tambien con el peso del vacio respecto al material estrucutral

% Materiales
ES.Es = 1; %Young del material estructural
ES.Ev = ES.gamma*ES.Es; % Young del material de los vacios
ES.nus=0.3; % Coeficiente de Poisson del material estructural.
ES.nuv=ES.nus; %Coeficiente de Poisson de los vacios

% Espesor
ES.esp=1; % Espesor de la estrucutra tridimensional equivalente

% Problema de elasticidad
% =1 si es Estado Plano de Tensiones
% =2 si es Estado Plano de Deformaciones
ES.PEL=1;

% Dimensiones y mallado
ES.Lx = 3.0;
ES.Ly = 1.0;
ES.Nx = 180; %Numero de subdivisiones segun X.
% NOTA: Para el de variaciones de tamaño de malla, no hacer /4 sino que
% hacer /3 y *3 (Eso por lo molesto del 0,55 del puente)
ES.Ny = ES.Nx/3; %Numero de subdivisiones segun Y

% Funcion de nivel inicial para los vacios
% Se define que si Psi>0 es material estructural
% y si es menor a 0 es material "vacio".
ES.Tolpsi = 1e-6; % tolerancia valores nodales de psi (no puede ser cero)
ES.fpsi = @(x,y) 1+0*x+0*y; % Todo material en primera instancia

ES.PropMat = [
    1   ES.Ev ES.nuv % Material vacio (siempre el 1)
    2   ES.Es ES.nus % Material estructural (siempre el 2)
    ];

% Matriz de Nodos
% [Nodo, X, Y]
xm = 0; xM = +ES.Lx;
ym = 0; yM = +ES.Ly;
dx = (xM-xm)/ES.Nx;
dy = (yM-ym)/ES.Ny;


ES.Nnodo = (ES.Nx+1)*(ES.Ny+1);
ES.Mnodo = zeros(ES.Nnodo,3);
pos = 0;
y = ym;
for iy = 1:ES.Ny+1
    x = xm;
    for ix = 1:ES.Nx+1
        pos = pos+1;
        ES.Mnodo(pos,:) = [pos, x, y];
        x = x + dx;
    end
    y = y + dy;
end

% Matriz conectividades elementos finitos
% IMPORTANTE INTENTAR METER SIMETRIA EN LA MALLA SI EL PROBLEMA ES
% SIMETRICO
% [Elem, Mat, Nod1, Nod2, Nod3]
ES.Nelem = 2*ES.Nx*ES.Ny;
ES.Melem = zeros(ES.Nelem,5);
pos = 0;
n1 = 1;
for iy = 1:ES.Ny
    for ix = 1:ES.Nx
        pos = pos+1;
        if ix<=ES.Nx/2
        ES.Melem(pos,1:5) = [pos, 2, n1 , n1+1, n1+2+ES.Nx];
        else
        ES.Melem(pos,1:5) = [pos, 2, n1,  n1+1, n1+1+ES.Nx];
        end
        pos = pos+1;
        if ix<=ES.Nx/2
        ES.Melem(pos,1:5) = [pos, 2, n1, n1+2+ES.Nx, n1+1+ES.Nx];
        else
        ES.Melem(pos,1:5) = [pos, 2, n1+1,  n1+2+ES.Nx, n1+1+ES.Nx];
        end
        n1 = n1 + 1;
    end
    % orden = not(orden);
    n1 = n1 + 1;
end

% Funcion phi inicial
ES.psi = ES.fpsi(ES.Mnodo(:,2), ES.Mnodo(:,3));

% =========================================================================
% == Condiciones de Contorno ==============================================
% =========================================================================

% NOTA:
% Todas las matrices que a continuacion se definen, aun en los casos
% en que no se coloque ninguno del tipo de apoyo que se seï¿½ala hay que
% indicar la matriz vacia.

% NOTA 2:
% Todo nodo que quede sin condicion de Dirichlet o Robin es considerado nodo
% con condicion de borde de Neumann. En el problema se resuelve el vector
% de fuerzas externas. Es mas que posible que haya nodos donde no se aclare
% el valor de la fuerza por ninguna fuente en cuyo caso esta es nula.

% Condiciones de Dirichlet - Apoyos Puntuales.
% ------------------------------------------------------
%
% Estos desplazamientos impuestos son solamente en los nodos.
% Pueden ser no nulos e inclinados.
%
% El "no aplicar" en bordes del elemento implica que no afectan a los
% grados de libertad extendidos.
%
% ATENCION: Si un nodo tiene cierto desplazamiento impuesto mas de una vez,
% el ultimo valor ingresado es el valido.
%
% Formato Matriz de nodos con todo el desplazamiento prescrito:
% ES.CB.Dir.Nod.Fijo = [Nodo, ValX, ValY]
% Nodo: Numero de nodo que posee desplazamiento prescrito en toda direccion
% ValX: Valor del desplazamiento impuesto segun x
% ValY: Valor del desplazamiento impuesto segun y
%
% Formato Matriz de nodos con apoyo deslizante.
% ES.CB.Dir.Nod.Desl = [Nodo, Angulo, Val];
% Nodo: Numero de nodo con desplazamiento prescrito en una direccion y
% libre en la ortogonal.
% Angulo: Direccion en la que el desplazamiento es libre. EN GRADOS.
%         Es positivo desde el eje x y antihorario.
%         NOTA: Angulos distintos de 0, 90, 180 o 270 ENLENTECEN MUCHO POR
%         TENER QUE ROTAR LA MATRIZ
% Val: Valor del desplazamiento prescrito en la direccion prescrita.
%      El sentido positivo de esta direccion es antihorario al eje definido
%      por angulo. Es decir si Angulo=20ï¿½, entonces el desplazamiento dado
%      por Val es positivo a un angulo de 110ï¿½.


ES.CB.Dir.Nod.Fijo=[[1;ES.Nx+1] , zeros(2,2) ];
ES.CB.Dir.Nod.Desl=[[(ES.Nx+1)*ES.Ny/2+1; (ES.Nx+1)*ES.Ny/2+ES.Nx+1 ] , zeros(2,2)]; %Idem

% Condiciones de Dirichlet - Apoyos En Bordes de Elemento
% -------------------------------------------------------
%
% %%%%%%%%%%%
% ATENCION!!!
%
% Imponer desplazamientos en bordes de elemento debe hacerse con cuidado.
% Si se imponen desplazamientos en dos bordes que comparten un nodo, el
% desplazamiento impuesto al nodo sera el del ultimo borde indicado en la
% lista de esta variable.
%
% EL USUARIO DEBE CUIDAR LA ADECUADA CONTINUIDAD!
% %%%%%%%%%%%
%
% Al apoyar el borde, puede verse afectado los grados de libertad
% extendidos. Se puede demostrar que si en el borde apoyado no hay
% interfaz, los grados de libertad extendidos siguen siendo libres, pero si
% hay interfaz los 2 extendidos, asociados al borde, son nulos (sin
% importar el valor de la imposicion)
%
% Para la numeracion de borde, se considera que el nodo inferior izquierdo
% es el primero y a partir de ahi se numeran los bordes en sentido
% antihorario.
%
% Dada la malla regular esto es:
%
%           2do
%     --------------              /|
%     |           /              / |
%     |          /              /  |
%     |         /              /   |
%     |        /              /    |
%     |       /              /     |
% 3ro |      / 1ro      3ro /      | 2do
%     |     /              /       |
%     |    /              /        |
%     |   /              /         |
%     |  /              /          |
%     | /              /           |
%     |/              /____________|
%                           1ro
%
% Formato Matriz de bordes con todo el desplazamiento prescrito:
% ES.CB.Dir.Lin.Fijo = [El, Bor, ValX, ValY]
% El: Elemento con el desplazamiento prescrito
% Bor: Numero de borde en el que se aplica
% ValX: Valor del desplazamiento impuesto segun x
% ValY: Valor del desplazamiento impuesto segun y
%
% Formato Matriz de bordes con apoyo deslizante.
% ES.CB.Dir.Lin.Desl = [El, Bor, Angulo, Val];
% El: Elemento con el desplazamiento prescrito en cierta direccion y libre
% en la ortogonal.
% Bor: Numero de borde en el que se aplica
% Angulo: Direccion en la que el desplazamiento es libre. EN GRADOS
%         NOTA: Angulos distintos de 0, 90, 180 o 270 ENLENTECEN MUCHO POR
%         TENER QUE ROTAR LA MATRIZ
% Val: Valor del desplazamiento prescrito en la direccion prescrita.
% Misma convencion de sentidos y angulos que para los nodales.

ES.CB.Dir.Lin.Fijo=[];

ES.CB.Dir.Lin.Desl=[]; %No tengo en este caso

% Condiciones de Neumann - Numero de Casos de carga
% -------------------------------------------------

ES.NLC=4;
% Si Este es mayor a 1 se deben generar las siguientes variables en formato
% cell, con tantas entradas como casos de carga.
% Si es =1 no pasa nada.
% Ademas hay que agregarle el "Cell" en el nombre

% Condiciones de Neumann - Fuerzas puntuales
%-------------------------------------------
%
% Formato:
% ES.CB.Neu.Punt = [Nodo, ValX, ValY]
% Nodo: Numero de nodo donde se aplica la fuerza
% ValX: Valor de la fuerza puntual segun X
% ValY: Valor de la fuerza puntual segun Y

ES.CB.NeuCell.Punt=cell(ES.NLC,1);

% Condiciones de Neumann - Fuerzas por unidad de longitud, uniformes, en bordes de elementos.
% -------------------------------------------------------------------------------------------
%
% Se seï¿½ala el numero de elemento y que borde es.
%
% Obs: Si se aplica una fuerza en un borde entre 2 elementos, se recomienda
% ingresarla en un solo elemento, en caso contrario la fuerza se duplica.
%
% La numeracion de los bordes es analoga a la ya descrita para el
% desplazamiento.
%
%
% Formato:
% ES.CB.Neu.Bord = [XoY, El, Bor, Val]
% XoY: Segun x (1) o Segun y (2)
% El: Elemento con la fuerza
% Bor: Numero de borde en la que se aplica
% Val: Valor de la fuerza por unidad de longitud

Linea=0.55/ES.Ly*ES.Ny;
if floor(Linea)~=Linea
    error('Revisa la cantidad de elementos segun Y porque la linea de carga no cae en borde de elemento')
end
Aux=ES.Nx*2*(Linea-1)+(2:2:2*ES.Nx);

ES.CB.NeuCell.Bord{2}=[2*ones(length(Aux),1), Aux', 2*ones(length(Aux),1),-10*ones(length(Aux),1)];
ES.CB.NeuCell.Bord{3}=[ones(length(Aux),1),Aux',2*ones(length(Aux),1),1*ones(length(Aux),1)];
ES.CB.NeuCell.Bord{4}=[ones(length(Aux),1),Aux',2*ones(length(Aux),1),-1*ones(length(Aux),1)];

% Condiciones de Neumann - Fuerzas de volumen (En general peso propio)
% --------------------------------------------------------------------
%
% Al asumirse que van a ser peso propio se va a considerar que TODOS los
% elementos la tienen y es uniforme en el dominio.
%
% Se admiten valores distintos para uno u otro material e incluso se puede
% imponer 0.
%
% Formato:
% ES.CB.Neu.Vol = [bx, by]
%
% bx: Componente segun el eje x de la fuerza de volumen
% by: Componente segun el eje y de la fuerza de volumen.
%
% SI O SI DEBE TENER DOS FILAS. LA PRIMERA SON LOS VACIOS Y LA SEGUNDA EL
% MATERIAL
ES.CB.NeuCell.Vol=cell(ES.NLC,1);

bMat=0.4;

ES.CB.NeuCell.Vol{1}=[0 -ES.gamma*bMat;
               0 -bMat];
for i=2:4
ES.CB.NeuCell.Vol{i}=zeros(2,2);
end
% Condiciones de Robin - Apoyos nodales elasticos
% -----------------------------------------------
%
% Estos apoyos elasticos se aplican solo en los nodos y pueden ser
% inclinados.
%
% Al ser solo en los nodos, no afectan a los grados de libertad extendidos
% pues estos generan desplazamiento nulo en el nodo.
%
% ATENCION: Si a un nodo se le definen mas de una condicion de apoyo
% elastico, los mismos se iran acumulando y esto probablemente sea
% incorrecto. EVITAR ESTO.
%
% Formato Matriz de nodos con toda direccion con apoyo elastico:
% ES.CB.Rob.Fijo = [Nodo, ValX, ValY]
% Nodo: Numero de nodo que posee apoyo elastico en toda direccion
% ValX: Valor de la constante elastica segun x (Unidades
% Fuerza/Desplazamiento)
% ValY: Valor de la constante elastica segun y (Mismas unidades)
%
% Formato Matriz de nodos con apoyo deslizante.
% ES.CB.Rob.Desl = [Nodo, Angulo, Val];
% Nodo: Numero de nodo con apoyo elastico en una direccion y
% libre en la ortogonal.
% Angulo: Direccion en la que el desplazamiento es libre. EN GRADOS.
%         Es positivo desde el eje x y antihorario.
%         NOTA: Angulos distintos de 0, 90, 180 o 270 ENLENTECEN MUCHO POR
%         TENER QUE ROTAR LA MATRIZ
% Val: Valor de la constante elastica en la direccion apoyada.

ES.CB.Rob.Fijo=[]; %No tengo de este tipo de apoyos en el ejemplo
ES.CB.Rob.Desl=[]; %Idem

% =============================================================================
% === Algoritmo de optimizacion ===============================================
% =============================================================================

% Variable logica que indica que puntos requieren tener level set estricto
% positivo: ( Material obligatorio )
ES.Must=ES.Mnodo(:,3)<=0.5501&ES.Mnodo(:,3)>=0.499;

ES.ChequeoSigma=ones(ES.Nelem,1); 
% Variable logica que indica en que elementos SI se requiere verificar las
% tensiones.

ES.MeshReg=true; % Es una malla regular (para seleccion del suavizado)


M = 0.3; % Porcentaje del volumen a buscar
c = 1200; %Valor para el coeficiente de penalidad del volumen al cuadrado (Lagrangiano fijo)
alpha = 0; %Valor para el coeficiente de penalidad de las tensiones.e
TolDeltaV = 0.005; % En fraccion, minimo deltavolumen si se usa el indicador
% como siguiente iteracion (Kappa=1)
% Si este volumen es muy chico se considera que se encontro la estructura
% optima y se modifica el lagrangiano.
DeltaVMax= 1; %Cuanto es el DeltaVolumen maximo admisible en una iteracion (En fraccion del total)
TolM = 1e-2; % Fraccion de volumen respecto a la cual M se considera que se llego
MaxSit=10;
Smooth=0.35;
DeCaSmooth = 1.2; % Coeficiente de decaimiento del smooth

% ESTA FORMA PUEDE NO SER LA MEJOR PARA MULTIPLES CASOS DE CARGA PUES
% PODRIA HABER UNO DE SUMA DONDE LA TENSION MAXIMA SUPERE EL LIMITE. SE
% PODRIA AJUSTAR CON LA TENSION DE FLUENCIA

ES = Optimizacion(ES,M,c,alpha,TolDeltaV,DeltaVMax,TolM,MaxSit,Smooth,DeCaSmooth);

save('EstructuraFinal.mat','ES','M','c','alpha','TolDeltaV','DeltaVMax','TolM','MaxSit','DeCaSmooth','Smooth','-v7.3')
