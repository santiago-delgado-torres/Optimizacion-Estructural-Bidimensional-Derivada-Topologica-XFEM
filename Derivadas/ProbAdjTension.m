function UAdj = ProbAdjTension(ES)
%
% Definicion y resolucion del problema adjunto al calculo de la derivada del
% funcional de tensiones segun Amstutz 2010.
%
% Devuelve UAdj, campo de "desplazamientos" del problema adjunto.

% =========================================================================
% === Deteccion de interfases =============================================
% =========================================================================

%ESAdj=interfase(ESAdj); % YA TA HECHO BOBO

% =========================================================================
% === Matriz de rigidez ===================================================
% =========================================================================

% La matriz de rigidez NO presenta cambios en el problema adjunto.
% ESAdj = matrizXFEM(ESAdj); YA ESTA HECHA BOBO


% =========================================================================
% === Armado de Vector de Fuerzas =========================================
% =========================================================================

% Rutina para hallar la matriz de rigidez tensionante
%
% =========================================================================
% === Inicializacion de variables para la matriz dispersa =================
% =========================================================================

KPos = 0; %Un indice de posicion que es comodo para ensamblar los vectores
TAM = 6 * 6 * ES.Nelem; % Si esto solo fuera FEM, cada elemento tiene 6 grados de libertad
% que lo afectan, entonces generaria una matriz de elemento 6 x 6.
% Entonces, pensando en matriz dispersa son 6 * 6 * la cantidad de
% elementos lo que hay que prealocar.
TAM = TAM + ( 18 * 18 - 6 * 6) * sum(ES.EI); %Pero ademas, cada elemento extendido tiene
% otros 12 grados de libertad que lo afectan. Entonces cada elemento
% extendido genera una matriz de elemento 18 x 18. 
% Se corrige entonces sacando (6*6*La Cantidad de Extendidos)
% Y sumando 18*18*La cantidad de extendidos.
% Ahora si esta correctamente definido el tamanio de la matriz de rigidez
% global

KVec = zeros(TAM,1); % Valores de las entradas
KLin = zeros(TAM,1); % Numeros de fila
KCol = zeros(TAM,1); % Numeros de columna

% =========================================================================
% === Iteracion en los elementos ==========================================
% =========================================================================


for ele = 1:ES.Nelem
    
    % -------------------------------------------
    % Variables que interesan a todo tipo -------
    
    % Nodos:
    ne = ES.Melem(ele,3:5);
    
    % Geometria del elemento:
    Xe = ES.Mnodo(ne,2);
    Ye = ES.Mnodo(ne,3); 
    
    % Derivadas de las funciones de forma segun coord. intrinsecas
    dN_detachi=[-1 1 0; -1 0 1];
    
    % Matriz Jacobiana y su determinante
    J = dN_detachi*[Xe,Ye];
    Jo = det(J);
    
    % Derivadas de las funciones de forma segun x e y.
    dN_dxy = J\dN_detachi; % Se hace con J\ en vez de lo usual de inv(J)* pues 
    % Matlab dice que asi es mas rapido.
    
    CharFun = ES.ChequeoSigma(ele); % Por ahora solo es DS

    
    
    
    % --------------------------------------------
    % Elemento extendido -------------------------
    if ES.EI(ele)
        
        ArX = ES.EGLX(ele,:); %Numero de aristas extendidas de este elemento
        
        if CharFun % Si no pertenece a DS ya ta, ni lo tomo en cuenta este y vale 0
            psie = ES.psi(ne); %Funcion de nivel para los nodos del elemento.

            ArX = ES.EGLX(ele,:); %Numero de aristas extendidas de este elemento
            NArX = ES.AriX(ArX,:); %Numero de nodo de esas aristas
            Uele = zeros(6,1); % Vector de desplazamientos nodales
            Uele(1:2:5) = ES.U(2*ne-1); % Se asignan los valores segun X
            Uele(2:2:6) = ES.U(2*ne); % Se asignan los valores segun Y        

            Uxele = zeros(12,1); % Vector de variables extendidas nodales
            ngx = full(ES.NodAriX(ne,ArX)); % Cuales son los grados de libertad segun x de los extendidos?
            Uxele(1:2:11) = ES.U(ngx); % Se asignan los valores segun x
            Uxele(2:2:12) = ES.U(ngx+1); % Y los segun y

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
            % segundo grado como es el caso)

            etaG = [0.5 0 0.5]; % Coordenadas eta de los puntos.
            xiG = [0.5 0.5 0]; % Coordenadas xi de los puntos
            PesoG = [1/6 1/6 1/6]; % Peso de los puntos.

            % INICIALIZACION DE MATRIZ DE RIGIDEZ DEL ELEMENTO

            Ke = zeros(18);

            % INTEGRACION EN SUBELEMENTOS
            %
            % Obs: Se dice subelementos como abuso de nomeclatura. En realidad
            % es siempre el mismo elemento y son subtriangulos para
            % integracion y definicion de funciones

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

                etase = etaSubElements(:,subEl);% Coordenadas eta de este subtriangulo.
                xise = xiSubElements(:,subEl); % Idem para la coordenada xi.
                Ar1 = Ar1SE(:,subEl); % Valores de Ar para g1
                Ar2 = Ar2SE(:,subEl); %Valores de Ar para g2

                JacSe = dN_detachi*[etase xise]; %Jacobiano del subelemento

                Jse = det( JacSe ); % Determinante jacobiano de este subtriangulo en el triangulo intrinseco

                dNse_detachi = JacSe\dN_detachi; % Derivada de las funciones de forma del subelemento
                % respecto a las coord. intrinsecas del elemento

                % Obtencion del material
                % Se obtiene cuanto vale el level set en un punto intermedio
                % del subelemento
                etaeval =  [1/3 1/3 1/3]* etase; % Interpolacion para el valor de eta.
                xieval = [1/3 1/3 1/3] *xise; %Idem para el valor de xi

                psiRep = [1-etaeval-xieval,etaeval,xieval]*psie; % Valor del level set representativo del SE

                % Seleccion del material

                if psiRep<0
                    E=ES.PropMat(1,2); %Young
                    nu=ES.PropMat(1,3); %Poisson
                    CharFunSE = 0; % Funcion caracteristica modificada si es vacio
                else
                    E=ES.PropMat(2,2); %Young
                    nu=ES.PropMat(2,3); %Poisson
                    CharFunSE = CharFun;
                end

                % Correccion si es EPD.

                if ES.PEL==1
                       mu = E/(2*(1+nu));
                       lambda = nu*E/(1-nu^2);
                    elseif ES.PEL==2
                       mu = E/(2*(1+nu));
                       lambda = nu*E/( ( 1+nu)*(1-2*nu));
                    end
                if ES.PEL==2
                    E = E/(1-nu^2);
                    nu = nu / (1-nu);
                end


                C = (E/(1-nu^2))* [ 1 nu 0;
                    nu 1 0;
                    0 0 0.5*(1-nu)]; % Matriz tensor elastico


                % Integracion en puntos de Gauss

                for puntoG=1:3 % Para los 3 puntos de Gauss anteriormente definidos

                    % Los puntos de Gauss son en sistemas de coordenadas
                    % intrinsecos a los subtriangulos.

                    % Para saber el valor de B se debe hayar los valores de eta
                    % y xi en el triangulo completo:
                    N1se = 1-etaG(puntoG)-xiG(puntoG) ; % Funcion de interpolacion clasica N1 en el subtriangulo
                    N2se = etaG(puntoG); %Idem funcion N2
                    N3se = xiG(puntoG); % Idem funcion N3
                    Nse = [N1se,N2se,N3se]; % Vector de estas 3 funciones
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
                    Sigma=C*Eps; 
                    
                    BTens = [ 4*mu + lambda ,   lambda - 2*mu  ,  0   ;
                          lambda - 2*mu ,   4*mu + lambda  ,  0   ;
                                0       ,        0         , 6*mu];
                            
                     BS = BTens*Sigma;
                            
                    SVM2 = 0.5 * ( BS(1)*Eps(1) + BS(2)*Eps(2) + BS(3)*Eps(3));

                    % Derivada de phi para el calculo de k1 y algo mas
                    t = SVM2 / (ES.SigmaM)^2; % Tension relativa
                    % Si t es tamaño "razonable" no hay problemas numéricos, pero si t pasa a
                    % ser muy grande la precisión es baja por lo cual se estudian los
                    % siguientes casos:
                    if t<= 1e7 % Se experimento y parecía funcionar hasta t = 1e9 pero para estar seguros
                        dPhidSVM2 =  (1 + t^32) ^ (-31/32) * t^31;
                    else % si se pasa, ya conviene asumir 1+t^32 = t^32 en cuyo caso
                        dPhidSVM2 = 1;
                    end % RECORDATORIO SE ELIMINO EL 1/SIGMAM2

                    % Valor de k1
                    k1 = dPhidSVM2 * CharFunSE / ES.SigmaM^2;
                    % Suma de este punto de Gauss

                    Ke = Ke + PesoG(puntoG) * Jo * Jse * ES.esp * k1 * B' * BTens* C * B ;

                end %En for en los puntos de Gauss


            end %Fin del for de recorrer los 4 subelementos.


            % Grados de libertad.

            DOF=zeros(18,1); %Inicializo

            DOF(1:2:5) = 2*ne'-1; % Grados de Libertad s/ x
            DOF(2:2:6) = 2*ne'; % Grados de libertad s/ y
            Aux = full(ES.NodAriX(ne,ArX)); %Se obtienen los grados de libertad segun x de los extendidos
            DOF(7:2:17) = Aux(:); % Grados de libertad extendidos s/ x
            DOF(8:2:18) = Aux(:)+1; % Grados de libertad extendidos s/ y 

            DOFMAT=repmat(DOF,1,18); %Se arma una matriz auxiliar de los grados de libertad para definir la matriz dispersa

            KLin( (KPos+1):(KPos+18^2) )=DOFMAT(:); % MATRIZ(:) recorre por columnas, entonces efectivamente esto va a devolver las filas donde se debe ubicar KVec.
            DOFMAT=transpose(DOFMAT); % Truco para hallar las columnas

            KCol( (KPos+1):(KPos+18^2) )=DOFMAT(:);

            KVec( (KPos+1):(KPos+18^2) )=Ke(:); 

            KPos=KPos+18^2;
        else
            % Grados de libertad.

            DOF=zeros(18,1); %Inicializo

            DOF(1:2:5) = 2*ne'-1; % Grados de Libertad s/ x
            DOF(2:2:6) = 2*ne'; % Grados de libertad s/ y
            Aux = full(ES.NodAriX(ne,ArX)); %Se obtienen los grados de libertad segun x de los extendidos
            DOF(7:2:17) = Aux(:); % Grados de libertad extendidos s/ x
            DOF(8:2:18) = Aux(:)+1; % Grados de libertad extendidos s/ y 

            DOFMAT=repmat(DOF,1,18); %Se arma una matriz auxiliar de los grados de libertad para definir la matriz dispersa

            KLin( (KPos+1):(KPos+18^2) )=DOFMAT(:); % MATRIZ(:) recorre por columnas, entonces efectivamente esto va a devolver las filas donde se debe ubicar KVec.
            DOFMAT=transpose(DOFMAT); % Truco para hallar las columnas

            KCol( (KPos+1):(KPos+18^2) )=DOFMAT(:);

            KVec( (KPos+1):(KPos+18^2) )=0; 

            KPos=KPos+18^2;
        end
    % --------------------------------------------
    % Elemento vacio -----------------------------
    elseif ES.VA(ele)
        
        Ke = 0; % El vacio tiene char fun nula
        
        % Grados de libertad.
        
        DOF=zeros(6,1); %Inicializo
        
        DOF(1:2:5) = 2*ne'-1; % Grados de Libertad s/ x
        DOF(2:2:6) = 2*ne'; % Grados de libertad s/ y
        
        DOFMAT=repmat(DOF,1,6); %Se arma una matriz auxiliar de los grados de libertad para definir la matriz dispersa

        KLin( (KPos+1):(KPos+36) )=DOFMAT(:); % MATRIZ(:) recorre por columnas, entonces efectivamente esto va a devolver las filas donde se debe ubicar KVec.
        DOFMAT=transpose(DOFMAT); % Truco para hallar las columnas

        KCol( (KPos+1):(KPos+36) )=DOFMAT(:);

        KVec( (KPos+1):(KPos+36) )=Ke(:); 

        KPos=KPos+36;
        
        
     % --------------------------------------------
     % Elemento tradicional -----------------------
    else 
        if CharFun
        % Asignacion de propiedades del material estructural
        E = ES.PropMat(2,2); %Young
        nu = ES.PropMat(2,3); %Poisson
         Uele = zeros(6,1); % Vector de desplazamientos nodales
        Uele(1:2:5) = ES.U(2*ne-1); % Se asignan los valores segun X
        Uele(2:2:6) = ES.U(2*ne); % Se asignan los valores segun Y
        % Correccion si es EPD.
        
        if ES.PEL==1
               mu = E/(2*(1+nu));
               lambda = nu*E/(1-nu^2);
            elseif ES.PEL==2
               mu = E/(2*(1+nu));
               lambda = nu*E/( ( 1+nu)*(1-2*nu));
            end
        if ES.PEL==2
            E = E/(1-nu^2);
            nu = nu / (1-nu);
        end

        
        C = (E/(1-nu^2))* [ 1 nu 0;
                        nu 1 0;
                        0 0 0.5*(1-nu)]; % Matriz tensor elastico
                
        % Matriz relacion desplazamiento-deformacion
        B = zeros(3,6);
        B(1,1:2:5)=dN_dxy(1,:);
        B(2,2:2:6)=dN_dxy(2,:);
        B(3,2:2:6)=dN_dxy(1,:);
        B(3,1:2:5)=dN_dxy(2,:);
        
        % Peso de Gauss. 
        % Nota: Se usa un unico punto de Gauss pues al ser B constante, 
        % y ser C, t, J tambien constantes, el integrando es constante y
        % por ende un unico punto integra bien.
        
        PesoG=0.5;
        Eps = B*Uele;
        Sigma=C*Eps;
        
        BTens = [ 4*mu + lambda ,   lambda - 2*mu  ,  0   ;
                          lambda - 2*mu ,   4*mu + lambda  ,  0   ;
                                0       ,        0         , 6*mu];
                            
         BS = BTens*Sigma;

        SVM2 = 0.5 * ( BS(1)*Eps(1) + BS(2)*Eps(2) + BS(3)*Eps(3));

        % Derivada de phi para el calculo de k1 y algo mas
        t = SVM2 / (ES.SigmaM)^2; % Tension relativa
        % Si t es tamaño "razonable" no hay problemas numéricos, pero si t pasa a
        % ser muy grande la precisión es baja por lo cual se estudian los
        % siguientes casos:
        if t<= 1e7 % Se experimento y parecía funcionar hasta t = 1e9 pero para estar seguros
            dPhidSVM2 =  (1 + t^32) ^ (-31/32) * t^31;
        else % si se pasa, ya conviene asumir 1+t^32 = t^32 en cuyo caso
            dPhidSVM2 = 1;
        end % RECORDATORIO SE ELIMINO EL 1/SIGMAM2

        % Valor de k1
        k1 = dPhidSVM2 / ES.SigmaM^2; % No preciso multiplicar por CharFun pues si era 0 ya no entrammos jeje
        % Matriz de rigidez
        Ke = PesoG * Jo * k1 * ES.esp * B' *BTens* C * B ;
        
        % Grados de libertad.
        
        DOF=zeros(6,1); %Inicializo
        
        DOF(1:2:5) = 2*ne'-1; % Grados de Libertad s/ x
        DOF(2:2:6) = 2*ne'; % Grados de libertad s/ y
        
        DOFMAT=repmat(DOF,1,6); %Se arma una matriz auxiliar de los grados de libertad para definir la matriz dispersa

        KLin( (KPos+1):(KPos+36) )=DOFMAT(:); % MATRIZ(:) recorre por columnas, entonces efectivamente esto va a devolver las filas donde se debe ubicar KVec.
        DOFMAT=transpose(DOFMAT); % Truco para hallar las columnas

        KCol( (KPos+1):(KPos+36) )=DOFMAT(:);

        KVec( (KPos+1):(KPos+36) )=Ke(:); 

        KPos=KPos+36;
        else
        DOF=zeros(6,1); %Inicializo
        
        DOF(1:2:5) = 2*ne'-1; % Grados de Libertad s/ x
        DOF(2:2:6) = 2*ne'; % Grados de libertad s/ y
        
        DOFMAT=repmat(DOF,1,6); %Se arma una matriz auxiliar de los grados de libertad para definir la matriz dispersa

        KLin( (KPos+1):(KPos+36) )=DOFMAT(:); % MATRIZ(:) recorre por columnas, entonces efectivamente esto va a devolver las filas donde se debe ubicar KVec.
        DOFMAT=transpose(DOFMAT); % Truco para hallar las columnas

        KCol( (KPos+1):(KPos+36) )=DOFMAT(:);

        KVec( (KPos+1):(KPos+36) )=0; 

        KPos=KPos+36;
        end
    end % End if en la seleccion de tipo de elemento
end %Fin de la iteracion en los elementos.


KAdj = sparse(KLin, KCol, KVec, 2*ES.Nnodo + ES.NGLX, 2*ES.Nnodo + ES.NGLX);

ES.Fext =  - KAdj * ES.U;

% =========================================================================
% === Prescripcion de Desplazamientos =====================================
% =========================================================================

% En el adjunto los desplazamientos son nulos no prescritos.
if ~isempty(ES.CB.Dir.Nod.Fijo)
ES.CB.Dir.Nod.Fijo(:,2:3)=0;
end
if ~isempty(ES.CB.Dir.Nod.Desl)
ES.CB.Dir.Nod.Desl(:,3)=0;
end
if ~isempty(ES.CB.Dir.Lin.Fijo)
ES.CB.Dir.Lin.Fijo(:,3:4)=0;
end
if ~isempty(ES.CB.Dir.Lin.Desl)
ES.CB.Dir.Lin.Desl(:,4)=0;
end

ES = DirichletXFEM(ES);

% =========================================================================
% === Resolucion del XFEM =================================================
% =========================================================================

% Reducir matrices

NoRest = 1:(2*ES.Nnodo + ES.NGLX ); %se inicializa vector de TODOS los grados de libertad
NoRest( ES.gdlfij==1 ) = []; % Pero se eliminan todas las entradas prescritas

Kred = ES.Ktotal( NoRest, NoRest ) ; %Se reduce la matriz de rigidez tangente
Fred = ES.Fext( NoRest ); %Idem para el vector de fuerzas

% Resolver sistea de ecuacioens

Ured = Kred\Fred; % Se resuelve las componentes no restringidas del vector de grados de libertad

% Reconstruir U

UAdj = zeros(2*ES.Nnodo + ES.NGLX , 1); %Se inicializa como nulo el vector de grados de libertad solucion
UAdj( ES.gdlfij==1 ) = ES.Upres( ES.gdlfij==1 ) ; % A los restrictos se les asigna su valor prescripto
UAdj(NoRest) = Ured; % y al resto el solucion del sistema de ecuaciones anterior

% Rotar aquellos que haya que rotar

for j = 1 : (ES.Nnodo + ES.NGLX/2) % Como se rotan de a, el numero de puntos a rotar  es la mitad de las filas
  
   theta = ES.GdLiRot(j); % Angulo que el nodo esta rotado
   if theta ~= 0
   Tr = [ cosd(theta) , -sind(theta) ;
             sind(theta) , cosd(theta) ]; %Matriz de rotacion

   UAdj( [j*2-1,j*2] ) = Tr*UAdj( [j*2-1,j*2] ); %Se devuelve a la base original
   end
end   


