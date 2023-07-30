function UAdj = ProbAdjTension(ESAdj)
%
% Definicion y resolucion del problema adjunto al calculo de la derivada del
% funcional de tensiones segun Amstutz 2010.
%
% Devuelve UAdj, campo de "desplazamientos" del problema adjunto.

% =========================================================================
% === Deteccion de interfases =============================================
% =========================================================================

ESAdj=interfase(ESAdj);

% =========================================================================
% === Matriz de rigidez ===================================================
% =========================================================================

% La matriz de rigidez NO presenta cambios en el problema adjunto.
ESAdj = matrizXFEM(ESAdj); 

% =========================================================================
% === Armado de Vector de Fuerzas =========================================
% =========================================================================

ESAdj.Fext = zeros(2*ESAdj.Nnodo + ESAdj.NGLX,1); % Se inicializan nuevamente como un vector nulo

E = ESAdj.PropMat(2,2); % Young
nu = ESAdj.PropMat(2,3); %nu
% Estos son los correctos ya que gamma0 del material es igual a 1; si no habria que
% dividir entre gamma0

% Parametros de Lame
if ESAdj.PEL==1
   mu = E/(2*(1+nu));
   lambda = nu*E/(1-nu^2);
elseif ESAdj.PEL==2
   mu = E/(2*(1+nu));
   lambda = nu*E/( ( 1+nu)*(1-2*nu));
end % Identico analisis para los de Lame

% Tensor B
BTens = [ 4*mu + lambda ,   lambda - 2*mu  ,  0   ;
      lambda - 2*mu ,   4*mu + lambda  ,  0   ;
            0       ,        0         , 6*mu];

% Tensor constitutivo
C = [lambda+2*mu, lambda,      0;
    lambda,       lambda+2*mu, 0;
    0,            0,           mu]; % Identico analisis; como sigma en el paper se calcula normalizado por gamma
% esto es lo mismo que calcular estos parametros solamente para el
% material.

% Llega la hora de generar el vector de fuerzas,
% mediante incroporar el efecto de las tensiones en todos los elementos
for ele = 1:ESAdj.Nelem
    % -------------------------------------------
    % Variables que interesan a todo tipo -------

    % Nodos:
    ne = ESAdj.Melem(ele,3:5);


    % Geometria del elemento:
    Xe = ESAdj.Mnodo(ne,2);
    Ye = ESAdj.Mnodo(ne,3);

    % Derivadas de las funciones de forma segun coord. intrinsecas
    dN_detachi=[-1 1 0; -1 0 1];

    % Matriz Jacobiana
    J = dN_detachi*[Xe,Ye];
    Jo = det(J);

    % Derivadas de las funciones de forma segun x e y.
    dN_dxy = J\dN_detachi; % Se hace con J\ en vez de lo usual de inv(J)* pues
    % Matlab dice que asi es mas rapido.
    
    FuncCaract = ESAdj.ChequeoSigma(ele); % Se obtiene la funcion caracteristica del dominio en este elemento
                

     % --------------------------------------------
    % Elemento extendido -------------------------
    if ESAdj.EI(ele)
        
        psie = ESAdj.psi(ne); %Funcion de nivel para los nodos del elemento.
        
        ArX = ESAdj.EGLX(ele,:); %Numero de aristas extendidas de este elemento
        NArX = ESAdj.AriX(ArX,:); %Numero de nodo de esas aristas
        
        Uele = zeros(6,1); % Vector de desplazamientos nodales
        Uele(1:2:5) = ESAdj.U(2*ne-1); % Se asignan los valores segun X
        Uele(2:2:6) = ESAdj.U(2*ne); % Se asignan los valores segun Y        

        Uxele = zeros(12,1); % Vector de variables extendidas nodales
        ngx = full(ESAdj.NodAriX(ne,ArX)); % Cuales son los grados de libertad segun x de los extendidos?
        Uxele(1:2:11) = ESAdj.U(ngx); % Se asignan los valores segun x
        Uxele(2:2:12) = ESAdj.U(ngx+1); % Y los segun y

        UGLele = [Uele;Uxele]; % Vector de todos los GL del elemento
        
        % PUNTOS INTERMEDIOS PARA DIVISION EN SUBELEMENTOS PARA INTEGRACION
        % Y DEFINICION DE LAS "g"

        % Punto entre nodo 1 y 2. En este xi = 0
        if psie(1)*psie(2)<0 
            etaR1 = psie(1)/(psie(1)-psie(2));

            if ( ne(1)==NArX(1,1) || ne(1)==NArX(1,2) ) && ( ne(2)==NArX(1,1) || ne(2) == NArX(1,2) ) % El punto pasa por la arista 1
                Ar1R1 = 1; % Indicador de que es el punto donde g1 vale 1
                Ar2R1 = 0; 
            else
                Ar2R1 = 1; %Inidcador de que es el punto donde g2 vale 1
                Ar1R1 = 0;
            end
        else
            etaR1 = 0.5;
            Ar1R1 = 0; Ar2R1 = 0;
        end

        % Punto entre nodo 2 y 3. En este xi = 1-eta
        if psie(2)*psie(3)<0
            etaR2 = psie(3)/(psie(3)-psie(2));
            if ( ne(2)==NArX(1,1) || ne(2)==NArX(1,2) ) && ( ne(3)==NArX(1,1) || ne(3) == NArX(1,2) ) % El punto pasa por la arista 1
                Ar1R2 = 1; % Indicador de que es el punto donde g1 vale 1
                Ar2R2 = 0; 
            else
                Ar2R2 = 1; %Inidcador de que es el punto donde g2 vale 1
                Ar1R2 = 0;
            end
        else
            etaR2 = 0.5;
            Ar1R2 = 0; Ar2R2 = 0;
        end  

        % Punto entre nodo 3 y 1. En este eta=0.        
        if psie(3)*psie(1)<0
            xiR3 = psie(1)/(psie(1)-psie(3));
            if ( ne(1)==NArX(1,1) || ne(1)==NArX(1,2) ) && ( ne(3)==NArX(1,1) || ne(3) == NArX(1,2) ) % El punto pasa por la arista 1
                Ar1R3 = 1; % Indicador de que es el punto donde g1 vale 1
                Ar2R3 = 0; 
            else
                Ar2R3 = 1; %Inidcador de que es el punto donde g2 vale 1
                Ar1R3 = 0;
            end
        else
            xiR3 = 0.5;
            Ar1R3 = 0; Ar2R3 = 0;
        end        
        
        % 3 PUNTOS DE GAUSS EN TRIANGULOS.
        % (3 puntos alcanza para integrar exactamente polinomios de
        % segundo grado como es aproximadamente el caso)

        etaG = [0.5 0 0.5]; % Coordenadas eta de los puntos.
        xiG = [0.5 0.5 0]; % Coordenadas xi de los puntos
        PesoG = [1/6 1/6 1/6]; % Peso de los puntos.

        % INTEGRACION EN SUBELEMENTOS
        %
        % Obs: Se dice subelementos como abuso de nomeclatura. En realidad
        % es siempre el mismo elemento y son subtriangulos para
        % integracion.

        etaSubElements = [0      1      0   etaR1;
            etaR1 etaR2   0   etaR2;
            0     etaR1 etaR2   0  ];
        % Ordenado en columnas la anterior variable es las coordenadas eta
        % de los subelementos.

        xiSubElements = [0      0       1      0    ;
            0   1-etaR2   xiR3  1-etaR2;
            xiR3   0    1-etaR2  xiR3  ];

        % Ordenado en columnas la anterior variable es las coordenadas xi
        % de los subelementos.    
        
        Ar1SE = [0      0      0     Ar1R1;
                 Ar1R1  Ar1R2  Ar1R3 Ar1R2;
                 Ar1R3  Ar1R1  Ar1R2 Ar1R3 ]; % Formato similar a los anteriores para los valores de g1 en los vertices de los subtriangulos
    
        Ar2SE = [0      0      0     Ar2R1;
                 Ar2R1  Ar2R2  Ar2R3 Ar2R2;
                 Ar2R3  Ar2R1  Ar2R2 Ar2R3 ]; % Formato similar a los anteriores para los valores de g2 en los vertices de los subtriangulos
             
             
        
        for subEl = 1:4 % Se recorren los 4 subelementos
            
            etase = etaSubElements(:,subEl); % Coordenadas eta de este subtriangulo.
            xise = xiSubElements(:,subEl); % Idem para la coordenada xi.
            Ar1 = Ar1SE(:,subEl); % Valores de Ar para g1
            Ar2 = Ar2SE(:,subEl); %Valores de Ar para g2
            
            JacSe = dN_detachi*[etase xise]; %Jacobiano del subelemento
            
            Jse = det( JacSe ); % Determinante jacobiano de este subtriangulo en el triangulo intrinseco
            
            dNse_detachi = JacSe\dN_detachi; % Derivada de las funciones de forma del subelemento
            

            % Obtencion del material
            % Se obtiene cuanto vale el level set en un punto intermedio
            % del subelemento
            etaeval =  [1/3 1/3 1/3]* etase; % Interpolacion para el valor de eta.
            xieval = [1/3 1/3 1/3] *xise; %Idem para el valor de xi
            
            psiRep = [1-etaeval-xieval,etaeval,xieval]*psie; % Valor del level set representativo del SE
            
            
            % Seleccion del material

            if psiRep<0 % Vacio
                gamma0 = ESAdj.gamma; %Valor de gamma
            else % Material
                gamma0 = 1; %Valor de gamma
            end

                                  
            % Integracion en puntos de Gauss

            for puntoG=1:3 % Para los 3 puntos de Gauss anteriormente definidos
                
                % Los puntos de Gauss son en sistemas de coordenadas
                % intrinsecos a los subtriangulos.

                % Para saber el valor de B se debe hallar los valores de eta
                % y xi en el triangulo completo:
                N1se = 1-etaG(puntoG) - xiG(puntoG); % Funcion de interpolacion clasica N1 en el subtriangulo
                N2se = etaG(puntoG); % Idem funcion N2
                N3se = xiG(puntoG); %Idem funcion N3
                Nse = [N1se,N2se,N3se]; % Vector de esas 3 funciones
                etaeval = Nse * etase; % Interpolacion para el valor de eta.
                xieval = Nse *xise; %Idem para el valor de xi
               
                % Matriz relacion desplazamiento-deformacion
                B = zeros(3,18);
                
                % Componentes clasicas
                B(1,1:2:5)=dN_dxy(1,:);
                B(2,2:2:6)=dN_dxy(2,:);
                B(3,2:2:6)=dN_dxy(1,:);
                B(3,1:2:5)=dN_dxy(2,:);
                                           
                % Para las componentes extendidas hay que hacer algun
                % calculo extra.
                
                % Funciones de interpolacion clasicas en este punto:
                N1 = 1-etaeval-xieval; % funcion de interpolacion clasica N1 para el elemento completo en este punto
                N2 = etaeval; %Idem N2
                N3 = xieval; %Idem N3
                N = [N1,N2,N3]; % vector de las mismas.
                
                
                % Funciones para definir los g
                g1 = Nse * Ar1;
                g2 = Nse * Ar2;
                
                % Funciones w
                w1 = N*g1;
                w2 = N*g2;
                
                % Derivadas de las funciones de interpolacion extendidas
                dg1_detachi = dNse_detachi*Ar1; % DErivada de g respecto a las intrinsecas
                dg2_detachi = dNse_detachi*Ar2;
                dw1_detachi = dg1_detachi*N + dN_detachi*g1; % Derivada de la funcion de interpolacion extendida respecto a eta y chi
                dw2_detachi = dg2_detachi*N + dN_detachi*g2;
                dw1_dxy = J\dw1_detachi; % Misma derivada pero respecto a x e y
                dw2_dxy = J\dw2_detachi;
                
                % Componentes extendidas
                B(1,7:2:11)=dw1_dxy(1,:);
                B(2,8:2:12)=dw1_dxy(2,:);
                B(3,8:2:12)=dw1_dxy(1,:);
                B(3,7:2:11)=dw1_dxy(2,:);
                
                B(1,13:2:17) = dw2_dxy(1,:);
                B(2,14:2:18) = dw2_dxy(2,:);
                B(3,14:2:18) = dw2_dxy(1,:);
                B(3,13:2:17) = dw2_dxy(2,:);


                % Vector de deformaciones elementales
                Eps = B *UGLele;
                
                % Vector de tensiones elementales
                Sigma=C*Eps; % Ya esta hecho la division entre gamma0 al dejar C fuera de la iteracion
                
                % Tension DE von Mises al cuadrado.
                SVM2 = ( Sigma(1) + Sigma(2) )^2 - 3*(Sigma(1)*Sigma(2) - Sigma(3)^2);

                % Derivada de phi para el calculo de k1
                SigmaM = ESAdj.SigmaM;
                t = SVM2 / (SigmaM)^2; % Tension relativa

                BS = BTens*Sigma; % Se aplica B del tensor
                
                % Calculo de k1:
                % Si t es tamaño "razonable" no hay problemas numéricos, pero si t pasa a
                % ser muy grande la precisión es baja por lo cual se estudian los
                % siguientes casos:
                if t<= 1e7 % Se experimento y parecía funcionar hasta t = 1e9 pero para estar seguros
                    k1 = FuncCaract*(1/SigmaM)^2 * (1 + t^32) ^ (-31/32) * t^31;
                else % si se pasa, ya conviene asumir 1+t^32 = t^32 en cuyo caso
                    k1 = FuncCaract * (1/SigmaM)^2;
                end            
 
                
                Felem = -PesoG(puntoG) * ESAdj.esp * Jo * Jse * B' * gamma0 * (k1 * BS); % Vector de fuerzas del elemento
                
                ESAdj.Fext( 2*ne - 1 ) = ESAdj.Fext( 2*ne - 1 ) + Felem(1:2:5); % Agregado a los Grados de libertad tradicionales s/ x
                ESAdj.Fext( 2*ne     ) = ESAdj.Fext( 2*ne     ) + Felem(2:2:6); % Agregado a los GL tradicionales s/ y
                ESAdj.Fext( full( ESAdj.NodAriX(ne,ArX(1)) ) ) = ESAdj.Fext( full( ESAdj.NodAriX(ne,ArX(1)) ) ) + Felem(7:2:11); % Agregado a los GL extendidos segun x
                ESAdj.Fext( full( ESAdj.NodAriX(ne,ArX(1)) )+1 ) = ESAdj.Fext( full( ESAdj.NodAriX(ne,ArX(1)) )+1 ) + Felem(8:2:12); % Agregado a los GL extendidos segun y
                ESAdj.Fext( full( ESAdj.NodAriX(ne,ArX(2)) ) ) = ESAdj.Fext( full( ESAdj.NodAriX(ne,ArX(2)) ) ) + Felem(13:2:17); % Agregado a los GL extendidos segun x
                ESAdj.Fext( full( ESAdj.NodAriX(ne,ArX(2)) )+1 ) = ESAdj.Fext( full( ESAdj.NodAriX(ne,ArX(2)) )+1 ) + Felem(14:2:18); % Agregado a los GL extendidos segun y

                
                
            end %En for en los puntos de Gauss
            
        end %Fin del for de recorrer los 4 subelementos.
        
        

     % --------------------------------------------
    % Elemento vacio -----------------------------
    elseif ESAdj.VA(ele)

             
        gamma0 = ESAdj.gamma; %En este elemento, si o si el valor de gamma0 es gamma.
        
        Uele = zeros(6,1); % Vector de desplazamientos nodales
        Uele(1:2:5) = ESAdj.U(2*ne-1); % Se asignan los valores segun X
        Uele(2:2:6) = ESAdj.U(2*ne); % Se asignan los valores segun Y
        
        

        % Matriz relacion desplazamiento-deformacion
        B = zeros(3,6);
        B(1,1:2:5)=dN_dxy(1,:);
        B(2,2:2:6)=dN_dxy(2,:);
        B(3,2:2:6)=dN_dxy(1,:);
        B(3,1:2:5)=dN_dxy(2,:);

        % Vector de deformaciones elementales
        Eps = B*Uele;
        
        % Vector de tensiones elementales
        Sigma=C*Eps; % Ya esta hecho la division entre gamma0 al dejar C fuera de la iteracion
        
        BS = BTens*Sigma; % Se aplica B del tensor
        
        % Tension de von Mises al cuadrado, calculada con la expresion de
        % BS.Epsilon
        SVM2 = 0.5 * ( BS(1)*Eps(1) + BS(2)*Eps(2) + BS(3)*Eps(3));

        % Derivada de phi para el calculo de k1
        SigmaM = ESAdj.SigmaM;
        t = SVM2 / (SigmaM)^2; % Tension relativa
        
        % Calculo de k1:
        % Si t es tamaño "razonable" no hay problemas numéricos, pero si t pasa a
        % ser muy grande la precisión es baja por lo cual se estudian los
        % siguientes casos:
        if t<= 1e7 % Se experimento y parecía funcionar hasta t = 1e9 pero para estar seguros
            k1 = FuncCaract*(1/SigmaM)^2 * (1 + t^32) ^ (-31/32) * t^31;
        else % si se pasa, ya conviene asumir 1+t^32 = t^32 en cuyo caso
            k1 = FuncCaract * (1/SigmaM)^2;
        end            

        
        % 1 PUNTO DE GAUSS EN TRIANGULOS.
        % (1 punto alcanza para integrar exactamente polinomios de
        % primer grado como es aproximadamente el caso)

        PesoG = 0.5; % Peso del punto.
        
        
        Felem = -PesoG * ESAdj.esp * Jo * B' * gamma0 * (k1 .* BS); % Vector de fuerzas del elemento
        
        % Incorporacion al vector de fuerzas
        
       
        ESAdj.Fext( 2*ne - 1 ) = ESAdj.Fext( 2*ne - 1 ) + Felem(1:2:5); % Agregado a los Grados de libertad tradicionales s/ x
        ESAdj.Fext( 2*ne     ) = ESAdj.Fext( 2*ne     ) + Felem(2:2:6); % Agregado a los GL tradicionales s/ y

        
        
     % --------------------------------------------
     % Elemento tradicional -----------------------
    else
                  
        
        Uele = zeros(6,1); % Vector de desplazamientos nodales
        Uele(1:2:5) = ESAdj.U(2*ne-1); % Se asignan los valores segun X
        Uele(2:2:6) = ESAdj.U(2*ne); % Se asignan los valores segun Y
        
        

        % Matriz relacion desplazamiento-deformacion
        B = zeros(3,6);
        B(1,1:2:5)=dN_dxy(1,:);
        B(2,2:2:6)=dN_dxy(2,:);
        B(3,2:2:6)=dN_dxy(1,:);
        B(3,1:2:5)=dN_dxy(2,:);

        % Vector de deformaciones elementales
        Eps = B*Uele;
        
        % Vector de tensiones elementales
        Sigma=C*Eps; % Recordar que ya tiene incorporada la normalizacion por gamma0 por la definicion de C.
        
        BS = BTens*Sigma; % Se aplica B del tensor
        
        % Tension de von Mises al cuadrado, calculada con la expresion de
        % BS.Epsilon
        SVM2 = 0.5 * ( BS(1)*Eps(1) + BS(2)*Eps(2) + BS(3)*Eps(3));

        % Derivada de phi para el calculo de k1
        SigmaM = ESAdj.SigmaM;
        t = SVM2 / (SigmaM)^2; % Tension relativa
        
        % Calculo de k1:
        % Si t es tamaño "razonable" no hay problemas numéricos, pero si t pasa a
        % ser muy grande la precisión es baja por lo cual se estudian los
        % siguientes casos:
        if t<= 1e7 % Se experimento y parecía funcionar hasta t = 1e9 pero para estar seguros
            k1 = FuncCaract*(1/SigmaM)^2 * (1 + t^32) ^ (-31/32) * t^31;
        else % si se pasa, ya conviene asumir 1+t^32 = t^32 en cuyo caso
            k1 = FuncCaract * (1/SigmaM)^2;
        end     
        
        % 1 PUNTO DE GAUSS EN TRIANGULOS.
        % (1 punto alcanza para integrar exactamente polinomios de
        % primer grado como es aproximadamente el caso)

        PesoG = 0.5; % Peso del punto.
        
        gamma = 1; % En este caso gamma = 1 pero no por eso deja de ser interesante considerarlo
        
        Felem =- PesoG * ESAdj.esp * Jo * B' * gamma * (k1 .* BS); % Vector de fuerzas del elemento
        
        % Incorporacion al vector de fuerzas
        
       
        ESAdj.Fext( 2*ne - 1 ) = ESAdj.Fext( 2*ne - 1 ) + Felem(1:2:5); % Agregado a los Grados de libertad tradicionales s/ x
        ESAdj.Fext( 2*ne     ) = ESAdj.Fext( 2*ne     ) + Felem(2:2:6); % Agregado a los GL tradicionales s/ y

        
        

    end % End if en la seleccion de tipo de elemento

end % Fin de recorrer todos los elementos

% =========================================================================
% === Agregado apoyos elasticos ===========================================
% =========================================================================

ESAdj = RobinXFEM(ESAdj); % Esto es una aproximacion!! Basada en que la rigidez se mantiene

% =========================================================================
% === Prescripcion de Desplazamientos =====================================
% =========================================================================

% En el adjunto los desplazamientos son nulos no prescritos.
if ~isempty(ESAdj.CB.Dir.Nod.Fijo)
ESAdj.CB.Dir.Nod.Fijo(:,2:3)=0;
end
if ~isempty(ESAdj.CB.Dir.Nod.Desl)
ESAdj.CB.Dir.Nod.Desl(:,3)=0;
end
if ~isempty(ESAdj.CB.Dir.Lin.Fijo)
ESAdj.CB.Dir.Lin.Fijo(:,3:4)=0;
end
if ~isempty(ESAdj.CB.Dir.Lin.Desl)
ESAdj.CB.Dir.Lin.Desl(:,4)=0;
end

ESAdj = DirichletXFEM(ESAdj);

% =========================================================================
% === Resolucion del XFEM =================================================
% =========================================================================

% Reducir matrices

NoRest = 1:(2*ESAdj.Nnodo + ESAdj.NGLX ); %se inicializa vector de TODOS los grados de libertad
NoRest( ESAdj.gdlfij==1 ) = []; % Pero se eliminan todas las entradas prescritas

Kred = ESAdj.Ktotal( NoRest, NoRest ) ; %Se reduce la matriz de rigidez tangente
Fred = ESAdj.Fext( NoRest ); %Idem para el vector de fuerzas

% Resolver sistea de ecuacioens

Ured = Kred\Fred; % Se resuelve las componentes no restringidas del vector de grados de libertad

% Reconstruir U

UAdj = zeros(2*ESAdj.Nnodo + ESAdj.NGLX , 1); %Se inicializa como nulo el vector de grados de libertad solucion
UAdj( ESAdj.gdlfij==1 ) = ESAdj.Upres( ESAdj.gdlfij==1 ) ; % A los restrictos se les asigna su valor prescripto
UAdj(NoRest) = Ured; % y al resto el solucion del sistema de ecuaciones anterior

% Rotar aquellos que haya que rotar

for j = 1 : (ESAdj.Nnodo + ESAdj.NGLX/2) % Como se rotan de a, el numero de puntos a rotar  es la mitad de las filas
  
   theta = ESAdj.GdLiRot(j); % Angulo que el nodo esta rotado
   if theta ~= 0
   Tr = [ cosd(theta) , -sind(theta) ;
             sind(theta) , cosd(theta) ]; %Matriz de rotacion

   UAdj( [j*2-1,j*2] ) = Tr*UAdj( [j*2-1,j*2] ); %Se devuelve a la base original
   end
end   


