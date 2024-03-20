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
ES.Problema = 'Portico';

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
ES.PEL=2;

% Dimensiones y mallado
ES.Lx = 20;
ES.luzX = 18;
ES.Ly = 9.0;
ES.luzY = 7.0;
ES.Nx = 15; %Numero de subdivisiones segun X por unidad de Lx
ES.Ny = ES.Nx; %Numero de subdivisiones segun Y por unidad de Ly

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
ES.xm = 0; ES.xM = +ES.Lx;
ES.ym = 0; ES.yM = +ES.Ly;
ES.dx = 1/ES.Nx;
ES.dy = 1/ES.Ny;



ES.Nnodo = (ES.Lx*ES.Nx+1)*(ES.Ly*ES.Ny+1);
ES.Mnodo = zeros(ES.Nnodo,3);
pos = 0;
y = ES.ym;
for iy = 1:ES.Ly*ES.Ny+1
    x = ES.xm;
    for ix = 1:ES.Lx*ES.Nx+1
        if y >= ES.luzY-1e-6 || abs(x-ES.Lx/2)>=ES.luzX/2-1e-6 % el 1e-66 es la clasica tolerancia
        pos = pos+1;
        ES.Mnodo(pos,:) = [pos, x, y];
        end
        x = x + ES.dx;
    end
    y = y + ES.dy;
end
ES.Mnodo=ES.Mnodo(1:(pos),:);
ES.Nnodo = length(ES.Mnodo(:,1));

% Matriz conectividades elementos finitos
% IMPORTANTE INTENTAR METER SIMETRIA EN LA MALLA SI EL PROBLEMA ES
% SIMETRICO
% [Elem, Mat, Nod1, Nod2, Nod3]
ES.Nelem = ((ES.Lx-ES.luzX)*ES.luzY+(ES.Lx)*(ES.Ly-ES.luzY))*ES.Nx*ES.Ny*2;
ES.Melem = zeros(ES.Nelem,5);
pos = 0;
n1 = 1;


% Columnas
for iy = 1:(ES.luzY*ES.Ny*2-1)
    if rem(iy,2)==1
    for ix = 1:ES.Nx*(ES.Lx-ES.luzX)/2
        pos = pos+1;
        ES.Melem(pos,1:5) = [pos, 2, n1 , n1+1, n1+1+((ES.Nx)*(ES.Lx-ES.luzX)/2 + 1)*2];
        pos = pos+1;
        ES.Melem(pos,1:5) = [pos, 2, n1, n1+1+((ES.Nx)*(ES.Lx-ES.luzX)/2 + 1)*2, n1+((ES.Nx)*(ES.Lx-ES.luzX)/2 + 1)*2];
        n1 = n1 + 1;
    end
    else
    for ix = 1:ES.Nx*(ES.Lx-ES.luzX)/2
        pos = pos+1;
        ES.Melem(pos,1:5) = [pos, 2, n1 , n1+1, n1+((ES.Nx)*(ES.Lx-ES.luzX)/2 + 1)*2];
        pos = pos+1;
        ES.Melem(pos,1:5) = [pos, 2, n1+1, n1+1+((ES.Nx)*(ES.Lx-ES.luzX)/2 + 1)*2, n1+((ES.Nx)*(ES.Lx-ES.luzX)/2 + 1)*2];
        n1 = n1 + 1;
    end
    end
        
    % orden = not(orden);
    n1 = n1 + 1;
end

for ix = 1:ES.Nx*(ES.Lx-ES.luzX)/2
        pos = pos+1;
        ES.Melem(pos,1:5) = [pos, 2, n1 , n1+1, n1+(ES.Nx)*(ES.Lx)+1];
        pos = pos+1;
        ES.Melem(pos,1:5) = [pos, 2, n1+1, n1+2+(ES.Nx)*(ES.Lx), n1+1+(ES.Nx)*(ES.Lx)];
        n1 = n1 + 1;
end

n1=n1+1;

% Tramo superior
for iy = 1:ES.Ny*(ES.Ly-ES.luzY)
    for ix = 1:ES.Nx*ES.Lx
        pos = pos+1;
        if ix<=ES.Nx*ES.Lx/2
        ES.Melem(pos,1:5) = [pos, 2, n1 , n1+1, n1+1+1+(ES.Nx)*(ES.Lx)];
        else
        ES.Melem(pos,1:5) = [pos, 2, n1,  n1+1, n1+1+(ES.Nx)*(ES.Lx)];
        end
        pos =pos+1;
        if ix<=ES.Nx*ES.Lx/2
        ES.Melem(pos,1:5) = [pos, 2, n1, n1+1+(ES.Nx)*(ES.Lx)+1, n1+(ES.Nx)*(ES.Lx)+1];
        else
        ES.Melem(pos,1:5) = [pos, 2, n1+1,  n1+1+(ES.Nx)*(ES.Lx)+1, n1+1+(ES.Nx)*(ES.Lx)];
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
% en que no se coloque ninguno del tipo de apoyo que se se�ala hay que
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
%      por angulo. Es decir si Angulo=20�, entonces el desplazamiento dado
%      por Val es positivo a un angulo de 110�.

ES.CB.Dir.Nod.Fijo=[]; %No tengo de este tipo de apoyos en el ejemplo
ES.CB.Dir.Nod.Desl=[]; %Idem

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

% Se apoya en toda direccion el borde de abajo a la izquierd ay a la
% derecha
Aux=(1):(2):ES.Nx*4; % Se puede ver que por la creacion de la malla son estos los elementos.

ES.CB.Dir.Lin.Fijo=[Aux',1*ones(length(Aux),1),zeros(length(Aux),1),zeros(length(Aux),1)];

ES.CB.Dir.Lin.Desl=[]; %No tengo en este caso

% Condiciones de Neumann - Numero de Casos de carga
% -------------------------------------------------

ES.NLC=3;
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

ES.CB.NeuCell.Punt=cell(3,1);

% Condiciones de Neumann - Fuerzas puntuales - A proyectar
%---------------------------------------------------------
%
% Formato:
% ES.CB.Neu.PuntProy = [Nodo, NodoDef, ValX, ValY]
% Nodo: Serie de numeros de nodos donde se podria aplicar la fuerza. Es un
% vector y se van recorriendo en orden.
% NodoDef: El nodo en el cual si ninuno de la lista sirvio... se aplica ahi
% jeje
% ValX: Valor de la fuerza puntual segun X
% ValY: Valor de la fuerza puntual segun Y
% Este si no lo quiero puedo ni crearlo

ES.CB.NeuCell.PuntProy=cell(3,1);

% Datos de Acero para usar de referencia jeje
% Pared que recibe (Cuando hagamos el combinado, hay que hacer los dos
% casos, viento de izquierda y viento de derecha.)
CWAceroBarlovento = 0.95;
CWAceroTecho = 0.80;
CWAceroSotavento = 0.70;
CMAcero = 1.7; 
SepPorticos= 4.5;
AnchoElementosY = 1/ES.Ny; % Esta parte va a ser algo mas compleja cuando introduzcamos mallas no regulares a los proticos jeje
AnchoElementosX = 1/ES.Nx;

CargaPuntualVientoBarlovento = CWAceroBarlovento*SepPorticos*AnchoElementosY;
CargaPuntualVientoSotavento = CWAceroSotavento*SepPorticos*AnchoElementosY;
CargaPuntualVientoTecho = CWAceroTecho * SepPorticos * AnchoElementosX;
CargaPuntualMuerta = CMAcero * SepPorticos * AnchoElementosX;

ES.CB.NeuCell.PuntProy{1} = cell(ES.luzX*ES.Nx-1,4);
ES.CB.NeuCell.PuntProy{2} = cell(ES.luzX*ES.Nx-1+2*(ES.luzY*ES.Ny),4);
ES.CB.NeuCell.PuntProy{3} = cell(ES.luzX*ES.Nx-1+2*(ES.luzY*ES.Ny),4);

% Pared izquierda

basico = 1;
 for nodis = 1:ES.luzY*ES.Ny
     Vec = zeros(1,(ES.Lx-ES.luzX)*ES.Nx/2+1);
     for i = 1:length(Vec)
         Vec(i) = (i-1) + (nodis-1)*2*((ES.Lx-ES.luzX)*ES.Nx/2+1 )+ basico;
     end
     Carg1 =  CargaPuntualVientoBarlovento;
     Carg2 =  - CargaPuntualVientoSotavento;
     ES.CB.NeuCell.PuntProy{2}{nodis,1}=Vec;
     ES.CB.NeuCell.PuntProy{3}{nodis,1}=Vec;
     ES.CB.NeuCell.PuntProy{2}{nodis,2}=Vec(1); % La aplico en la cara izquierda si no tengo donde
     ES.CB.NeuCell.PuntProy{3}{nodis,2} = Vec(1);
     ES.CB.NeuCell.PuntProy{2}{nodis,3}=Carg1;
     ES.CB.NeuCell.PuntProy{3}{nodis,3}=Carg2;
     ES.CB.NeuCell.PuntProy{2}{nodis,4}=0;
     ES.CB.NeuCell.PuntProy{3}{nodis,4}=0;
 end
 

% Techo

basico1 = find( abs(ES.Mnodo(:,2)-1)<1e-6 & abs(ES.Mnodo(:,3)-9) <1e-6);
basico2 = find( abs(ES.Mnodo(:,2)-1)<1e-6 & abs(ES.Mnodo(:,3)-7) <1e-6);

 for nodis = 1:ES.luzX*ES.Nx-1
     Vec1 = zeros(1,(ES.Ly-ES.luzY)*ES.Ny+1); % LA de viento la aplico arriba y la de la carga muerta abajo
     Vec2 = zeros(1,(ES.Ly-ES.luzY)*ES.Ny+1);
     for i = 1:length(Vec1)
         Vec1(i) = -(ES.Lx*ES.Nx+1)* (i-1) + nodis + basico1;
         Vec2(i) = (ES.Lx*ES.Nx+1)* (i-1) + nodis + basico2;
     end
     Def1 = nodis + basico1; % La carga de viento la voy a aplicar arriba si no tengo donde
     Def2 = nodis + basico2; % La carga muerta la voy a aplicar abajo si no tengo donde
     Carg1 =  CargaPuntualVientoTecho;
     Carg2 = - CargaPuntualMuerta;
     ES.CB.NeuCell.PuntProy{2}{nodis+ES.luzY*ES.Ny,1}=Vec1;
     ES.CB.NeuCell.PuntProy{3}{nodis+ES.luzY*ES.Ny,1}=Vec1;
     ES.CB.NeuCell.PuntProy{1}{nodis,1}=Vec2;
     
     ES.CB.NeuCell.PuntProy{2}{nodis+ES.luzY*ES.Ny,2}=Def1;
     ES.CB.NeuCell.PuntProy{3}{nodis+ES.luzY*ES.Ny,2}=Def1;
     ES.CB.NeuCell.PuntProy{1}{nodis,2}=Def2;
     
     ES.CB.NeuCell.PuntProy{2}{nodis+ES.luzY*ES.Ny,3}=0;
     ES.CB.NeuCell.PuntProy{3}{nodis+ES.luzY*ES.Ny,3}=0;
     ES.CB.NeuCell.PuntProy{1}{nodis,3}=0;
     
     ES.CB.NeuCell.PuntProy{2}{nodis+ES.luzY*ES.Ny,4}=Carg1;
     ES.CB.NeuCell.PuntProy{3}{nodis+ES.luzY*ES.Ny,4}=Carg1;
     ES.CB.NeuCell.PuntProy{1}{nodis,4}=Carg2;
 end
 

basico = 1;

 for nodis = 1:ES.luzY*ES.Ny
     Vec =zeros(1,(ES.Lx-ES.luzX)*ES.Nx/2+1);
     for i = 1:length(Vec)
        Vec(i) =  (nodis)*2*((ES.Lx-ES.luzX)*ES.Nx/2+1 )+ basico - i;
     end
     Carg1 =  CargaPuntualVientoSotavento;
     Carg2 = -CargaPuntualVientoBarlovento;
     ES.CB.NeuCell.PuntProy{2}{nodis+ES.luzY*ES.Ny+ES.luzX*ES.Nx-1,1}=Vec;
     ES.CB.NeuCell.PuntProy{3}{nodis+ES.luzY*ES.Ny+ES.luzX*ES.Nx-1,1}=Vec;
     ES.CB.NeuCell.PuntProy{2}{nodis+ES.luzY*ES.Ny+ES.luzX*ES.Nx-1,2}=Vec(1); % La aplico en la cara derecha si no tengo donde
     ES.CB.NeuCell.PuntProy{3}{nodis+ES.luzY*ES.Ny+ES.luzX*ES.Nx-1,2} = Vec(1);
     ES.CB.NeuCell.PuntProy{2}{nodis+ES.luzY*ES.Ny+ES.luzX*ES.Nx-1,3}=Carg1;
     ES.CB.NeuCell.PuntProy{3}{nodis+ES.luzY*ES.Ny+ES.luzX*ES.Nx-1,3}=Carg2;
     ES.CB.NeuCell.PuntProy{2}{nodis+ES.luzY*ES.Ny+ES.luzX*ES.Nx-1,4}=0;
     ES.CB.NeuCell.PuntProy{3}{nodis+ES.luzY*ES.Ny+ES.luzX*ES.Nx-1,4}=0;
 end
% Condiciones de Neumann - Fuerzas por unidad de longitud, uniformes, en bordes de elementos.
% -------------------------------------------------------------------------------------------
%
% Se se�ala el numero de elemento y que borde es.
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


ES.CB.NeuCell.Bord=cell(3,1);

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

bMat=0;
ES.CB.NeuCell.Vol=cell(3,1);
ES.CB.NeuCell.Vol{1}=[0 -ES.gamma*bMat;
    0 -bMat];
ES.CB.NeuCell.Vol{2}=zeros(2,2);
ES.CB.NeuCell.Vol{3}=zeros(2,2);

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
ES.Must=zeros(ES.Nnodo,1); 

ES.ChequeoSigma=ones(ES.Nelem,1); 
% Variable logica que indica en que elementos SI se requiere verificar las
% tensiones.

ES.MeshReg=true; % Es una malla regular (para seleccion del suavizado)


M = 0.7; % Porcentaje del volumen a buscar
c = 300; %Valor para el coeficiente de penalidad del volumen al cuadrado (Lagrangiano fijo)
alpha = 0; %Valor para el coeficiente de penalidad de las tensiones.e
TolDeltaV = 0.005; % En fraccion, minimo deltavolumen si se usa el indicador
% como siguiente iteracion (Kappa=1)
% Si este volumen es muy chico se considera que se encontro la estructura
% optima y se modifica el lagrangiano.
DeltaVMax= 1; %Cuanto es el DeltaVolumen maximo admisible en una iteracion (En fraccion del total)
TolM = 1e-2; % Fraccion de volumen respecto a la cual M se considera que se llego
MaxSit=10;
Smooth=0.4;
DeCaSmooth = 1.2; % Coeficiente de decaimiento del smooth

% ESTA FORMA PUEDE NO SER LA MEJOR PARA MULTIPLES CASOS DE CARGA PUES
% PODRIA HABER UNO DE SUMA DONDE LA TENSION MAXIMA SUPERE EL LIMITE. SE
% PODRIA AJUSTAR CON LA TENSION DE FLUENCIA

ES = Optimizacion(ES,M,c,alpha,TolDeltaV,DeltaVMax,TolM,MaxSit,Smooth,DeCaSmooth);

save('EstructuraFinal.mat','ES','M','c','alpha','TolDeltaV','DeltaVMax','TolM','MaxSit','DeCaSmooth','Smooth','-v7.3')

