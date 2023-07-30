% =========================================================================
% == Inicio limpio de terminal ============================================
% =========================================================================
clc
clear all
close all

% =========================================================================
% == Definicion de la Estructura ==========================================
% =========================================================================
%
% A seguir se define la estructura, se deben utilizar UNIDADES COMPATIBLES.
% NO DEBEN modificarse los nombres de los campos de la variable ES.

% Nombre del problema
ES.Problema = 'ProblemaEjemplo';

% Materiales
ES.Ev = 1e-6; % Young del material de los vacios
ES.Es = 1; %Young del material estructural 
ES.nuv=0; %Coeficiente de Poisson de los vacios
ES.nus=0; % Coeficiente de Poisson del material estructural.

% Espesor
ES.esp=1; % Espesor de la estrucutra tridimensional equivalente

% Problema de elasticidad
% =1 si es Estado Plano de Tensiones
% =2 si es Estado Plano de Deformaciones
ES.PEL=1;

% Dimensiones y mallado
ES.Lx = 1.0;
ES.Ly = 1;
ES.Nx = 50; %Numero de subdivisiones segun X. (Cantidad de elementos)
ES.Ny = ES.Nx; %Numero de subdivisiones segun Y

% Funcion de nivel inicial para los vacios 
% Se define que si Psi>0 es material estructural
% y si es menor a 0 es material "vacio".
ES.Tolpsi = 1e-6; % tolerancia valores nodales de psi (no puede ser cero)
ES.fpsi = @(x,y) (sqrt((x-0).^2 + (y-0).^2) - 0.23); % Bola
%ES.fpsi=@(x,y) 0.49/2-x;

ES.PropMat = [
    1   ES.Ev ES.nuv % Material vacio (siempre el 1)
    2   ES.Es ES.nus % Material estructural (siempre el 2)
    ];

% Matriz de Nodos
% [Nodo, X, Y]
xm =-ES.Lx/2; xM = +ES.Lx/2;
ym = -ES.Ly/2; yM = +ES.Ly/2;
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
% [Elem, Mat, Nod1, Nod2, Nod3]
ES.Nelem = 2*ES.Nx*ES.Ny;
ES.Melem = zeros(ES.Nelem,5);
pos = 0;
n1 = 1;
for iy = 1:ES.Ny
    for ix = 1:ES.Nx
        pos = pos+1;
        ES.Melem(pos,1:5) = [pos, 2, n1 , n1+1, n1+2+ES.Nx];
        pos = pos+1;
        ES.Melem(pos,1:5) = [pos, 2, n1, n1+2+ES.Nx, n1+1+ES.Nx];
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

S = 1;
R = 0.5;
ep = .23;
al = ES.Ev/ES.Es;
g = al;
E = ES.Ev;
v = ES.nuv;
muu = E/2/(1+v);
la = v*E/(1-v^2);

AR0 = (R^2*g*la + 2*R^2*g*muu)/(4*(R^2*muu - ep^2*muu + R^2*g*la + R^2*g*muu + ep^2*g*muu));
AR1 = -(R^4*ep^2*g^2*la^2 + 3*R^4*ep^2*g^2*la*muu + 2*R^4*ep^2*g^2*muu^2 - R^4*ep^2*g*la^2 - 3*R^4*ep^2*g*la*muu - 2*R^4*ep^2*g*muu^2 - R^2*ep^4*g^2*la^2 - 3*R^2*ep^4*g^2*la*muu - 2*R^2*ep^4*g^2*muu^2 + R^2*ep^4*g*la^2 + 3*R^2*ep^4*g*la*muu + 2*R^2*ep^4*g*muu^2)/(R^8*g^2*la^2 + 4*R^8*g^2*la*muu + 3*R^8*g^2*muu^2 + 2*R^8*g*la^2 + 8*R^8*g*la*muu + 10*R^8*g*muu^2 + R^8*la^2 + 4*R^8*la*muu + 3*R^8*muu^2 + 4*R^6*ep^2*g^2*la^2 + 8*R^6*ep^2*g^2*la*muu + 4*R^6*ep^2*g^2*muu^2 + 8*R^6*ep^2*g*la*muu + 8*R^6*ep^2*g*muu^2 - 4*R^6*ep^2*la^2 - 16*R^6*ep^2*la*muu - 12*R^6*ep^2*muu^2 - 6*R^4*ep^4*g^2*la^2 - 12*R^4*ep^4*g^2*la*muu - 6*R^4*ep^4*g^2*muu^2 - 12*R^4*ep^4*g*la*muu - 12*R^4*ep^4*g*muu^2 + 6*R^4*ep^4*la^2 + 24*R^4*ep^4*la*muu + 18*R^4*ep^4*muu^2 + 4*R^2*ep^6*g^2*la^2 + 12*R^2*ep^6*g^2*la*muu + 12*R^2*ep^6*g^2*muu^2 + 4*R^2*ep^6*g*la*muu - 4*R^2*ep^6*la^2 - 16*R^2*ep^6*la*muu - 12*R^2*ep^6*muu^2 + ep^8*g^2*la^2 + 4*ep^8*g^2*la*muu + 3*ep^8*g^2*muu^2 - 2*ep^8*g*la^2 - 8*ep^8*g*la*muu - 6*ep^8*g*muu^2 + ep^8*la^2 + 4*ep^8*la*muu + 3*ep^8*muu^2);
AR2 = (R^2*muu + R^2*g*la + R^2*g*muu)/(4*(R^2*muu - ep^2*muu + R^2*g*la + R^2*g*muu + ep^2*g*muu));
AR3 = -(R^4*ep^2*g^2*la^2 + 2*R^4*ep^2*g^2*la*muu + R^4*ep^2*g^2*muu^2 + 2*R^4*ep^2*g*la*muu + 2*R^4*ep^2*g*muu^2 - R^4*ep^2*la^2 - 4*R^4*ep^2*la*muu - 3*R^4*ep^2*muu^2 - R^2*ep^4*g^2*la^2 - 2*R^2*ep^4*g^2*la*muu - R^2*ep^4*g^2*muu^2 - 2*R^2*ep^4*g*la*muu - 2*R^2*ep^4*g*muu^2 + R^2*ep^4*la^2 + 4*R^2*ep^4*la*muu + 3*R^2*ep^4*muu^2)/(2*(R^8*g^2*la^2 + 4*R^8*g^2*la*muu + 3*R^8*g^2*muu^2 + 2*R^8*g*la^2 + 8*R^8*g*la*muu + 10*R^8*g*muu^2 + R^8*la^2 + 4*R^8*la*muu + 3*R^8*muu^2 + 4*R^6*ep^2*g^2*la^2 + 8*R^6*ep^2*g^2*la*muu + 4*R^6*ep^2*g^2*muu^2 + 8*R^6*ep^2*g*la*muu + 8*R^6*ep^2*g*muu^2 - 4*R^6*ep^2*la^2 - 16*R^6*ep^2*la*muu - 12*R^6*ep^2*muu^2 - 6*R^4*ep^4*g^2*la^2 - 12*R^4*ep^4*g^2*la*muu - 6*R^4*ep^4*g^2*muu^2 - 12*R^4*ep^4*g*la*muu - 12*R^4*ep^4*g*muu^2 + 6*R^4*ep^4*la^2 + 24*R^4*ep^4*la*muu + 18*R^4*ep^4*muu^2 + 4*R^2*ep^6*g^2*la^2 + 12*R^2*ep^6*g^2*la*muu + 12*R^2*ep^6*g^2*muu^2 + 4*R^2*ep^6*g*la*muu - 4*R^2*ep^6*la^2 - 16*R^2*ep^6*la*muu - 12*R^2*ep^6*muu^2 + ep^8*g^2*la^2 + 4*ep^8*g^2*la*muu + 3*ep^8*g^2*muu^2 - 2*ep^8*g*la^2 - 8*ep^8*g*la*muu - 6*ep^8*g*muu^2 + ep^8*la^2 + 4*ep^8*la*muu + 3*ep^8*muu^2));
AR4 = -(ep^2*(R^8*g^2*la^2 + 2*R^8*g^2*la*muu + R^8*g^2*muu^2 + 2*R^8*g*la*muu + 2*R^8*g*muu^2 - R^8*la^2 - 4*R^8*la*muu - 3*R^8*muu^2 + R^2*ep^6*g^2*la^2 + 4*R^2*ep^6*g^2*la*muu + 3*R^2*ep^6*g^2*muu^2 - 2*R^2*ep^6*g*la^2 - 8*R^2*ep^6*g*la*muu - 6*R^2*ep^6*g*muu^2 + R^2*ep^6*la^2 + 4*R^2*ep^6*la*muu + 3*R^2*ep^6*muu^2))/(2*(R^8*g^2*la^2 + 4*R^8*g^2*la*muu + 3*R^8*g^2*muu^2 + 2*R^8*g*la^2 + 8*R^8*g*la*muu + 10*R^8*g*muu^2 + R^8*la^2 + 4*R^8*la*muu + 3*R^8*muu^2 + 4*R^6*ep^2*g^2*la^2 + 8*R^6*ep^2*g^2*la*muu + 4*R^6*ep^2*g^2*muu^2 + 8*R^6*ep^2*g*la*muu + 8*R^6*ep^2*g*muu^2 - 4*R^6*ep^2*la^2 - 16*R^6*ep^2*la*muu - 12*R^6*ep^2*muu^2 - 6*R^4*ep^4*g^2*la^2 - 12*R^4*ep^4*g^2*la*muu - 6*R^4*ep^4*g^2*muu^2 - 12*R^4*ep^4*g*la*muu - 12*R^4*ep^4*g*muu^2 + 6*R^4*ep^4*la^2 + 24*R^4*ep^4*la*muu + 18*R^4*ep^4*muu^2 + 4*R^2*ep^6*g^2*la^2 + 12*R^2*ep^6*g^2*la*muu + 12*R^2*ep^6*g^2*muu^2 + 4*R^2*ep^6*g*la*muu - 4*R^2*ep^6*la^2 - 16*R^2*ep^6*la*muu - 12*R^2*ep^6*muu^2 + ep^8*g^2*la^2 + 4*ep^8*g^2*la*muu + 3*ep^8*g^2*muu^2 - 2*ep^8*g*la^2 - 8*ep^8*g*la*muu - 6*ep^8*g*muu^2 + ep^8*la^2 + 4*ep^8*la*muu + 3*ep^8*muu^2));
BR0 = -(R^2*(R^6*g^2*la^2 + 3*R^6*g^2*la*muu + 2*R^6*g^2*muu^2 + R^6*g*la^2 + 5*R^6*g*la*muu + 6*R^6*g*muu^2 - 3*R^2*ep^4*g^2*la^2 - 9*R^2*ep^4*g^2*la*muu - 6*R^2*ep^4*g^2*muu^2 + 3*R^2*ep^4*g*la^2 + 9*R^2*ep^4*g*la*muu + 6*R^2*ep^4*g*muu^2 + 4*ep^6*g^2*la^2 + 14*ep^6*g^2*la*muu + 12*ep^6*g^2*muu^2 - 4*ep^6*g*la^2 - 14*ep^6*g*la*muu - 12*ep^6*g*muu^2))/(R^8*g^2*la^2 + 4*R^8*g^2*la*muu + 3*R^8*g^2*muu^2 + 2*R^8*g*la^2 + 8*R^8*g*la*muu + 10*R^8*g*muu^2 + R^8*la^2 + 4*R^8*la*muu + 3*R^8*muu^2 + 4*R^6*ep^2*g^2*la^2 + 8*R^6*ep^2*g^2*la*muu + 4*R^6*ep^2*g^2*muu^2 + 8*R^6*ep^2*g*la*muu + 8*R^6*ep^2*g*muu^2 - 4*R^6*ep^2*la^2 - 16*R^6*ep^2*la*muu - 12*R^6*ep^2*muu^2 - 6*R^4*ep^4*g^2*la^2 - 12*R^4*ep^4*g^2*la*muu - 6*R^4*ep^4*g^2*muu^2 - 12*R^4*ep^4*g*la*muu - 12*R^4*ep^4*g*muu^2 + 6*R^4*ep^4*la^2 + 24*R^4*ep^4*la*muu + 18*R^4*ep^4*muu^2 + 4*R^2*ep^6*g^2*la^2 + 12*R^2*ep^6*g^2*la*muu + 12*R^2*ep^6*g^2*muu^2 + 4*R^2*ep^6*g*la*muu - 4*R^2*ep^6*la^2 - 16*R^2*ep^6*la*muu - 12*R^2*ep^6*muu^2 + ep^8*g^2*la^2 + 4*ep^8*g^2*la*muu + 3*ep^8*g^2*muu^2 - 2*ep^8*g*la^2 - 8*ep^8*g*la*muu - 6*ep^8*g*muu^2 + ep^8*la^2 + 4*ep^8*la*muu + 3*ep^8*muu^2);
BR2 = -(R^2*(R^6*g^2*la^2 + 4*R^6*g^2*la*muu + 3*R^6*g^2*muu^2 + 2*R^6*g*la^2 + 8*R^6*g*la*muu + 10*R^6*g*muu^2 + R^6*la^2 + 4*R^6*la*muu + 3*R^6*muu^2 - 3*R^2*ep^4*g^2*la^2 - 6*R^2*ep^4*g^2*la*muu - 3*R^2*ep^4*g^2*muu^2 - 6*R^2*ep^4*g*la*muu - 6*R^2*ep^4*g*muu^2 + 3*R^2*ep^4*la^2 + 12*R^2*ep^4*la*muu + 9*R^2*ep^4*muu^2 + 4*ep^6*g^2*la^2 + 12*ep^6*g^2*la*muu + 12*ep^6*g^2*muu^2 + 4*ep^6*g*la*muu - 4*ep^6*la^2 - 16*ep^6*la*muu - 12*ep^6*muu^2))/(2*(R^8*g^2*la^2 + 4*R^8*g^2*la*muu + 3*R^8*g^2*muu^2 + 2*R^8*g*la^2 + 8*R^8*g*la*muu + 10*R^8*g*muu^2 + R^8*la^2 + 4*R^8*la*muu + 3*R^8*muu^2 + 4*R^6*ep^2*g^2*la^2 + 8*R^6*ep^2*g^2*la*muu + 4*R^6*ep^2*g^2*muu^2 + 8*R^6*ep^2*g*la*muu + 8*R^6*ep^2*g*muu^2 - 4*R^6*ep^2*la^2 - 16*R^6*ep^2*la*muu - 12*R^6*ep^2*muu^2 - 6*R^4*ep^4*g^2*la^2 - 12*R^4*ep^4*g^2*la*muu - 6*R^4*ep^4*g^2*muu^2 - 12*R^4*ep^4*g*la*muu - 12*R^4*ep^4*g*muu^2 + 6*R^4*ep^4*la^2 + 24*R^4*ep^4*la*muu + 18*R^4*ep^4*muu^2 + 4*R^2*ep^6*g^2*la^2 + 12*R^2*ep^6*g^2*la*muu + 12*R^2*ep^6*g^2*muu^2 + 4*R^2*ep^6*g*la*muu - 4*R^2*ep^6*la^2 - 16*R^2*ep^6*la*muu - 12*R^2*ep^6*muu^2 + ep^8*g^2*la^2 + 4*ep^8*g^2*la*muu + 3*ep^8*g^2*muu^2 - 2*ep^8*g*la^2 - 8*ep^8*g*la*muu - 6*ep^8*g*muu^2 + ep^8*la^2 + 4*ep^8*la*muu + 3*ep^8*muu^2));
BR3 = -(R^2*(ep^2*muu - ep^2*g*muu))/(2*(R^2*muu - ep^2*muu + R^2*g*la + R^2*g*muu + ep^2*g*muu));
BR4 = -(R^4*(g - 1)*(ep^8*g*la^2 - 3*ep^8*muu^2 - ep^8*la^2 + 3*ep^8*g*muu^2 + R^4*ep^4*la^2 + 3*R^4*ep^4*muu^2 - 4*ep^8*la*muu + R^4*ep^4*g*la^2 + R^4*ep^4*g*muu^2 + 4*ep^8*g*la*muu + 4*R^4*ep^4*la*muu + 2*R^4*ep^4*g*la*muu))/(2*(R^8*g^2*la^2 + 4*R^8*g^2*la*muu + 3*R^8*g^2*muu^2 + 2*R^8*g*la^2 + 8*R^8*g*la*muu + 10*R^8*g*muu^2 + R^8*la^2 + 4*R^8*la*muu + 3*R^8*muu^2 + 4*R^6*ep^2*g^2*la^2 + 8*R^6*ep^2*g^2*la*muu + 4*R^6*ep^2*g^2*muu^2 + 8*R^6*ep^2*g*la*muu + 8*R^6*ep^2*g*muu^2 - 4*R^6*ep^2*la^2 - 16*R^6*ep^2*la*muu - 12*R^6*ep^2*muu^2 - 6*R^4*ep^4*g^2*la^2 - 12*R^4*ep^4*g^2*la*muu - 6*R^4*ep^4*g^2*muu^2 - 12*R^4*ep^4*g*la*muu - 12*R^4*ep^4*g*muu^2 + 6*R^4*ep^4*la^2 + 24*R^4*ep^4*la*muu + 18*R^4*ep^4*muu^2 + 4*R^2*ep^6*g^2*la^2 + 12*R^2*ep^6*g^2*la*muu + 12*R^2*ep^6*g^2*muu^2 + 4*R^2*ep^6*g*la*muu - 4*R^2*ep^6*la^2 - 16*R^2*ep^6*la*muu - 12*R^2*ep^6*muu^2 + ep^8*g^2*la^2 + 4*ep^8*g^2*la*muu + 3*ep^8*g^2*muu^2 - 2*ep^8*g*la^2 - 8*ep^8*g*la*muu - 6*ep^8*g*muu^2 + ep^8*la^2 + 4*ep^8*la*muu + 3*ep^8*muu^2));
 


Aux=[1:1:(ES.Nx+1),2*(ES.Nx+1):(ES.Nx+1):ES.Nnodo,(ES.Nnodo-1):-1:((ES.Nx+1)*ES.Ny+1),((ES.Nx+1)*(ES.Ny-1)+1):-(ES.Nx+1):(ES.Nx+2)];
XImp=0*Aux;
YImp=0*Aux;
for i=1:length(Aux)
 xeste = ES.Mnodo(Aux(i),2);
 yeste = ES.Mnodo(Aux(i),3);
 r=sqrt(xeste.^2+yeste^2);
 t = atan2(yeste,xeste); if t<0; t=t+2*pi; end
    
    ur2 = -al*(S*(BR4*la*cos(2*t) + BR4*muu*cos(2*t) + BR3*la*r.^2 - 2*AR2*muu*r.^4 + BR3*muu*r.^2 - 2*AR4*la*r.^2.*cos(2*t) + 2*AR3*la*r.^6.*cos(2*t) + BR2*la*r.^4.*cos(2*t) - 4*AR4*muu*r.^2.*cos(2*t) + BR2*muu*r.^4.*cos(2*t)))./(2*muu*r.^3*(la + muu));
ut2 = al*(S*sin(2*t).*(4*AR3*la*r.^6 - BR4*muu - BR4*la + BR2*la*r.^4 - 2*AR4*muu*r.^2 + 6*AR3*muu*r.^6 + BR2*muu*r.^4))./(2*muu*r.^3*(la + muu));

ux2 = cos(t).*ur2 - sin(t).*ut2;
uy2 = sin(t).*ur2 + cos(t).*ut2;
    XImp(i)=ux2;
    YImp(i)=uy2;
end

ES.CB.Dir.Nod.Fijo=[Aux',XImp',YImp']; %No tengo de este tipo de apoyos en el ejemplo
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

ES.CB.Dir.Lin.Fijo=[];

ES.CB.Dir.Lin.Desl=[];

% Condiciones de Neumann - Fuerzas puntuales 
%-------------------------------------------
%
% Formato:
% ES.CB.Neu.Punt = [Nodo, ValX, ValY]
% Nodo: Numero de nodo donde se aplica la fuerza
% ValX: Valor de la fuerza puntual segun X
% ValY: Valor de la fuerza puntual segun Y

ES.CB.Neu.Punt=[]; %En este caso no vamos a colocar ninguna puntual :)

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

% En este ejemplo se aplica una fuerza distribuida de valor q
% en todo el borde derecho y hacia la derecha.

ES.CB.Neu.Bord=[];

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

ES.CB.Neu.Vol=[0 0;
               0 0]; 
           
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



% =========================================================================
% = Solucion del problema  XFEM ===========================================
% =========================================================================
ES = xfem(ES);

% =========================================================================
% = Figuras ===============================================================
% =========================================================================

% Mapa de UX
figure
figXFEM(ES,'UX',1,1)

% Mapa de UY
figure
figXFEM(ES,'UY',1,1)

% Mapa de Deformada
figure
figXFEM(ES,'DEF',0,0,0.3)



% Mapa de tensiones segun x
figure 
figXFEM(ES,'SigmaX',1,1)

% Mapa de tensiones segun y
figure
figXFEM(ES,'SigmaY',1,1)

% MApa de tensiones segun XY
figure
figXFEM(ES,'SigmaXY',1,0)

% Mapa de tensiones de VonMises
figure
figXFEM(ES,'SigmaMises',1,1)

% Solo el borde
figure 
figXFEM(ES,'B'); %  Esta es interesante acoplarla con otras.

% Mapa de Psi
figure
figXFEM(ES,'Psi',1,1)
% Nota: En el de Psi no importa el valor de Vacum. Siempre se van a dibujar los vacios porque si no no tiene sentido

% Un plot negro material blanco vacio
figure
figXFEM(ES,'Mat',1,1) 
% Tampoco importa el valor de vacum.












