function ES = DerTop_Complacencia_XFEM(ES)

% Rutina para hallar la derivada topologica de la complacencia
% segun el paper de Lopes 2015. Recibe un x2 respecto a la del paper pues
% la del paper es La energia potencial total que es -0.5 Complacencia.
%
% Calcuando las derivadas especiales en los elementos finitos y luego
% asignando al nodo un valor ponderado (VP) segun la siguiente formula:
%
% VP = 6 * \int N * DT.
% Siendo N la funcion de forma del nodo en el elemento, DT la derivada
% topologica y la integral en el elemento de referencia ( (0,0) (1,0) (0,1) )
% El valor 6 es porque si DT fuera cte = 1 quisieramos que VP = 1 y la
% integral de N *1 = 1/6
%
% Como algunos (la mayoria) de los nodos tienen mas de un elemento, se
% calcula el VP de cada elemento. Luego para definir la DT del
% nodo se calcula la media de los VP.
%
% Devuelve o actualiza:
%
% ES.DTC:     Derivada topologica de la complacencia. Tiene tanta entradas como nodo la malla

% =========================================================================
% === Primeras variables ==================================================
% =========================================================================

SumVP = zeros(ES.Nnodo,1); % Se inicializa un vector de la suma de los VP.
EleCounter = zeros(ES.Nnodo,1); % Y un vector que cuenta los elementos por nodo  para calcular la media.


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

    % Matriz Jacobiana
    J = dN_detachi*[Xe,Ye];

    % Derivadas de las funciones de forma segun x e y.
    dN_dxy = J\dN_detachi; % Se hace con J\ en vez de lo usual de inv(J)* pues
    % Matlab dice que asi es mas rapido.

     % --------------------------------------------
    % Elemento extendido -------------------------
    if ES.EI(ele)
       
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

        % 4 PUNTOS DE GAUSS EN TRIANGULOS.
        % (4 puntos alcanza para integrar exactamente polinomios de
        % tercer grado como es el caso)

        etaG = [1/3 0.6 0.2 0.2]; % Coordenadas eta de los puntos.
        xiG = [1/3 0.2 0.6 0.2]; % Coordenadas xi de los puntos
        PesoG = [-0.28125 0.260416666666666667 0.260416666666666667  0.260416666666666667]; % Peso de los puntos.

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
            % respecto a las coord. intrinsecas del elemento
            
            % Obtencion del material
            % Se obtiene cuanto vale el level set en un punto intermedio
            % del subelemento
            etaeval =  [1/3 1/3 1/3]* etase; % Interpolacion para el valor de eta.
            xieval = [1/3 1/3 1/3] *xise; %Idem para el valor de xi
            
            psiRep = [1-etaeval-xieval,etaeval,xieval]*psie; % Valor del level set representativo del SE
            
            if psiRep<0
                E=ES.PropMat(1,2); %Young
                nu=ES.PropMat(1,3); %Poisson
                gammaP = 1/ES.gamma ; %Para el material vacio el tensor de polarizacion es si le asigno un material de 1/gamma veces la rigidez.
                b=ES.CB.Neu.Vol(1,:);  % Fuerza de volumen
            else
                E=ES.PropMat(2,2); %Young
                nu=ES.PropMat(2,3); %Poisson
                gammaP = ES.gamma ; %Para el material estructural el tensor de polarizacion es si le asigno un material de gamma veces la rigidez.
                b=ES.CB.Neu.Vol(2,:);  % Fuerza de volumen
            end

            if ES.PEL==1
               mu = E/(2*(1+nu));
               lambda = nu*E/(1-nu^2);
            elseif ES.PEL==2
               mu = E/(2*(1+nu));
               lambda = nu*E/( ( 1+nu)*(1-2*nu));
               E = E/(1-nu^2);
               nu = nu / (1-nu);
            end

            C = (E/(1-nu^2))* [ 1 nu 0;
                            nu 1 0;
                            0 0 0.5*(1-nu)]; % Matriz tensor elastico

            % Alpha 1 y 2 de la expresion de Lopes
            alpha1= (lambda+mu)/mu;
            alpha2= (lambda+3*mu)/(lambda+mu);
            

            % Ya aprovechemos y sumemos al contador:
            EleCounter(ne) = EleCounter(ne) + Jse.*ones(3,1);

            % Integracion en puntos de Gauss

            for puntoG=1:4 % Para los 4 puntos de Gauss anteriormente definidos

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
                 %En realidad la componente del tensor fuera de la diagonal es la
                 %mitad de la asignada en ese vector a la tercer componente
                % pero como luego se va a hacer un producto escalar con un tensor simetrico
                % podemos hacer este producto una vez usando este valor.
                % Ademas C se arma considerando que se esta considerando la
                % distorsion angular y no la componente del tensor.

                 % Vector de tensiones elementales
                Sigma=C*Eps;

                % Se pasa a formato matricial (es mas comodo para el tensor de
                % polarizacion)
                Sigma=[Sigma(1), Sigma(3); Sigma(3), Sigma(2)];

                trS = Sigma(1,1) + Sigma(2,2); %Traza del tensor

                 % Formato matricial del tensor de polarizacion aplicado a Sigma
                SigmaP = 0.5*( (1-gammaP)/(1+gammaP*alpha2) ) * ( (1+alpha2)*Sigma + 0.5*(alpha1-alpha2)*( (1-gammaP)/(1+gammaP*alpha1) )*trS*eye(2) );

                % Pasaje a formato vectorial
                SigmaP = [SigmaP(1,1);SigmaP(2,2);SigmaP(1,2)];

                Ux = [ N , w1 , w2 ] * UGLele(1:2:17) ; % Componente x del desplazamiento interpolado al punto de Gauss
                Uy =  [ N , w1 , w2 ] * UGLele(2:2:18) ; % Idem componetne Y

                DerTopG = SigmaP'*Eps + (1-gammaP) * (b(1) * Ux + b(2)*Uy ); % Derivada topologica en el punto de Gauss
                % Nota, como b ya tiene incorporado el gamma por ser el del
                % vacio y el gammaP ya tiene el 1/gamma . esta formula es
                % correcta
                DerTopG = 2 * DerTopG; %Pasaje de ptencial total a complacencia.

                ValPerNodo =  ( N * DerTopG )'; % Valores por nodo en este punto

                SumVP(ne) = SumVP(ne)+ 6 * Jse *PesoG(puntoG) * ValPerNodo; % Se suma este punto de Gauss



            end %En for en los puntos de Gauss


        end %Fin del for de recorrer los 4 subelementos.


     % --------------------------------------------
    % Elemento vacio -----------------------------
    elseif ES.VA(ele)

        % Ya aprovechemos y sumemos al contador:
        EleCounter(ne) = EleCounter(ne) + ones(3,1);

        % Asignacion de propiedades del material vacio
        E = ES.PropMat(1,2); %Young
        nu = ES.PropMat(1,3); %Poisson

        if ES.PEL==1
           mu = E/(2*(1+nu));
           lambda = nu*E/(1-nu^2);
        elseif ES.PEL==2
           mu = E/(2*(1+nu));
           lambda = nu*E/( ( 1+nu)*(1-2*nu));
           E = E/(1-nu^2);
           nu = nu / (1-nu);
        end

        Uele = zeros(6,1); % Vector de desplazamientos nodales
        Uele(1:2:5) = ES.U(2*ne-1); % Se asignan los valores segun X
        Uele(2:2:6) = ES.U(2*ne); % Se asignan los valores segun Y

        % Matriz relacion desplazamiento-deformacion
        B = zeros(3,6);
        B(1,1:2:5)=dN_dxy(1,:);
        B(2,2:2:6)=dN_dxy(2,:);
        B(3,2:2:6)=dN_dxy(1,:);
        B(3,1:2:5)=dN_dxy(2,:);

        % Vector de deformaciones elementales
        Eps = B*Uele;
         %En realidad la componente del tensor fuera de la diagonal es la
         %mitad de la asignada en ese vector a la tercer componente
        % pero como luego se va a hacer un producto escalar con un tensor simetrico
        % podemos hacer este producto una vez usando este valor.
        % Ademas C se arma considerando que se esta considerando la
        % distorsion angular y no la componente del tensor.


        C = (E/(1-nu^2))* [ 1 nu 0;
                        nu 1 0;
                        0 0 0.5*(1-nu)]; % Matriz tensor elastico

        % Vector de tensiones elementales
        Sigma=C*Eps;

        % Se pasa a formato matricial (es mas comodo para el tensor de
        % polarizacion)
        Sigma=[Sigma(1), Sigma(3); Sigma(3), Sigma(2)];

        trS = Sigma(1,1) + Sigma(2,2); %Traza del tensor

        gammaP = 1/ES.gamma ; %Para el material vacio el tensor de polarizacion es si le asigno un material de 1/gamma veces la rigidez.

        % Alpha 1 y 2 de la expresion de Lopes
        alpha1= (lambda+mu)/mu;
        alpha2= (lambda+3*mu)/(lambda+mu);

        % Formato matricial del tensor de polarizacion aplicado a Sigma
        SigmaP = 0.5*( (1-gammaP)/(1+gammaP*alpha2) ) * ( (1+alpha2)*Sigma + 0.5*(alpha1-alpha2)*( (1-gammaP)/(1+gammaP*alpha1) )*trS*eye(2) );

        % Pasaje a formato vectorial
        SigmaP = [SigmaP(1,1);SigmaP(2,2);SigmaP(1,2)];

        % Fuerza de volumen
        b=ES.CB.Neu.Vol(1,:);


        % 3 PUNTOS DE GAUSS EN TRIANGULOS.
        % (3 puntos alcanza para integrar exactamente polinomios de
        % segundo grado como es el caso)

        etaG = [0.5 0 0.5]; % Coordenadas eta de los puntos.
        xiG = [0.5 0.5 0]; % Coordenadas xi de los puntos
        PesoG = [1/6 1/6 1/6]; % Peso de los puntos.

        % Integracion en puntos de Gauss

        for puntoG = 1:3 %Para los 3 puntos de Gauss anteriormente definidos

            N1 = 1-etaG(puntoG) -xiG(puntoG) ; % Funcion de forma N1 clasica para el elemento en este punto
            N2 = etaG(puntoG); %Idem N2
            N3 = xiG(puntoG); %Idem N3
            N=[N1,N2,N3]; % En formato vector estas 3

            Ux = N * Uele(1:2:5) ; % Componente x del desplazamiento interpolado al punto de Gauss
            Uy = N * Uele(2:2:6) ; % Idem componetne Y

            DerTopG = dot(SigmaP, Eps) + (1-gammaP) * (b(1) * Ux + b(2)*Uy ); % Derivada topologica en el punto de Gauss
            % Nota, como b ya tiene incorporado el gamma por ser el del
            % vacio y el gammaP ya tiene el 1/gamma . esta formula es
            % correcta
            DerTopG = 2 * DerTopG; %Pasaje de ptencial total a complacencia.

            ValPerNodo = ( N * DerTopG )'; % Valores por nodo en este punto

            SumVP(ne) =SumVP(ne)+ 6 * PesoG(puntoG) * ValPerNodo; % Se suma este punto de Gauss

        end % Fin de integrar en puntos de gauss

      % --------------------------------------------
     % Elemento tradicional -----------------------
    else

        % Ya aprovechemos y sumemos al contador:
        EleCounter(ne) = EleCounter(ne) + ones(3,1);

        % Asignacion de propiedades del material estructural
        E = ES.PropMat(2,2); %Young
        nu = ES.PropMat(2,3); %Poisson

        if ES.PEL==1
           mu = E/(2*(1+nu));
           lambda = nu*E/(1-nu^2);
        elseif ES.PEL==2
           mu = E/(2*(1+nu));
           lambda = nu*E/( ( 1+nu)*(1-2*nu));
           E = E/(1-nu^2);
           nu = nu / (1-nu);
        end

        Uele = zeros(6,1); % Vector de desplazamientos nodales
        Uele(1:2:5) = ES.U(2*ne-1); % Se asignan los valores segun X
        Uele(2:2:6) = ES.U(2*ne); % Se asignan los valores segun Y

        % Matriz relacion desplazamiento-deformacion
        B = zeros(3,6);
        B(1,1:2:5)=dN_dxy(1,:);
        B(2,2:2:6)=dN_dxy(2,:);
        B(3,2:2:6)=dN_dxy(1,:);
        B(3,1:2:5)=dN_dxy(2,:);

        % Vector de deformaciones elementales
        Eps = B*Uele;
         %En realidad la componente del tensor fuera de la diagonal es la
         %mitad de la asignada en ese vector a la tercer componente
        % pero como luego se va a hacer un producto escalar con un tensor simetrico
        % podemos hacer este producto una vez usando este valor.
        % Ademas C se arma considerando que se esta considerando la
        % distorsion angular y no la componente del tensor.


        C = (E/(1-nu^2))* [ 1 nu 0;
                        nu 1 0;
                        0 0 0.5*(1-nu)]; % Matriz tensor elastico

        % Vector de tensiones elementales
        Sigma=C*Eps;

        % Se pasa a formato matricial (es mas comodo para el tensor de
        % polarizacion)
        Sigma=[Sigma(1), Sigma(3); Sigma(3), Sigma(2)];

        trS = Sigma(1,1) + Sigma(2,2); %Traza del tensor

        gammaP = ES.gamma ; %Para el material estructural el tensor de polarizacion es si le asigno un material de gamma veces la rigidez.

        % Alpha 1 y 2 de la expresion de Lopes
        alpha1= (lambda+mu)/mu;
        alpha2= (lambda+3*mu)/(lambda+mu);

        % Formato matricial del tensor de polarizacion aplicado a Sigma
        SigmaP = 0.5*( (1-gammaP)/(1+gammaP*alpha2) ) * ( (1+alpha2)*Sigma + 0.5*(alpha1-alpha2)*( (1-gammaP)/(1+gammaP*alpha1) )*trS*eye(2) );

        % Pasaje a formato vectorial
        SigmaP = [SigmaP(1,1);SigmaP(2,2);SigmaP(1,2)];

        % Fuerza de volumen
        b=ES.CB.Neu.Vol(2,:);


        % 3 PUNTOS DE GAUSS EN TRIANGULOS.
        % (3 puntos alcanza para integrar exactamente polinomios de
        % segundo grado como es el caso)

        etaG = [0.5 0 0.5]; % Coordenadas eta de los puntos.
        xiG = [0.5 0.5 0]; % Coordenadas xi de los puntos
        PesoG = [1/6 1/6 1/6]; % Peso de los puntos.

        % Integracion en puntos de Gauss

        for puntoG = 1:3 %Para los 3 puntos de Gauss anteriormente definidos

            N1 = 1-etaG(puntoG) -xiG(puntoG) ; % Funcion de forma N1 clasica para el elemento en este punto
            N2 = etaG(puntoG); %Idem N2
            N3 = xiG(puntoG); %Idem N3
            N=[N1,N2,N3]; % En formato vector estas 3

            Ux = N  * Uele(1:2:5) ; % Componente x del desplazamiento interpolado al punto de Gauss
            Uy = N  * Uele(2:2:6) ; % Idem componetne Y

            DerTopG = dot(SigmaP, Eps) + (1-gammaP) * (b(1) * Ux + b(2)*Uy ); % Derivada topologica en el punto de Gauss
            DerTopG = 2 * DerTopG; %Pasaje de ptencial total a complacencia.
            
            ValPerNodo = ( N * DerTopG )'; % Valores por nodo en este punto

            SumVP(ne) = SumVP(ne)+6 * PesoG(puntoG) * ValPerNodo; % Se suma este punto de Gauss

        end % Fin de integrar en puntos de gauss

    end % End if en la seleccion de tipo de elemento

end % Fin de recorrer todos los elementos


ES.DTC = SumVP./EleCounter; % Se calcula la media
