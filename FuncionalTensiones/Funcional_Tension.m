function FuncTension = Funcional_Tension(ES)
% Rutina para hallar el funcional de la limitacion de tensiones
% segun el paper de Amstutz 2010. 
%
% Se calcula el valor del integrando dentro de los elementos finitos y se
% va sumando mediante integración en puntos de Gauss.
%
% Observación: Para el valor usado de p (ver paper) el integrando es casi
% lineal con las tensiones de von mises al cuadrado.
% Como las tensiones son uniformes en el elemento tradicional también lo es
% el integrando en estos. Se puede calcular con un UNICO punto de Gauss.
% En el caso del elemento extendido las tensiones son lineales, por lo
% tanto el cuadrado de von mises es cuadratico por lo que PHI es aprox
% cuadratico lo que requiere 3 puntos de Gauss.
%
%
% Devuelve:
%
% FuncTension:  El valor del funcional de tensiones para el estado actual.
%

% =========================================================================
% === Primeras variables ==================================================
% =========================================================================

FuncTension = 0;

% =========================================================================
% === Iteracion en los elementos ==========================================
% =========================================================================

for ele = 1:ES.Nelem
    if ES.ChequeoSigma(ele) % Si no se debe chequear sigma ni hagamos las cuentas para este elemento
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
                % respecto a las coord. intrinsecas del elemento
                
                
                % Obtencion del material
                % Se obtiene cuanto vale el level set en un punto intermedio
                % del subelemento
                etaeval =  [1/3 1/3 1/3]* etase; % Interpolacion para el valor de eta.
                xieval = [1/3 1/3 1/3] *xise; %Idem para el valor de xi

                psiRep = [1-etaeval-xieval,etaeval,xieval]*psie; % Valor del level set representativo del SE

                % Seleccion del material
                
                E=ES.PropMat(2,2); %Young 
                nu=ES.PropMat(2,3); %Poisson % Es material estructural siempre segun la definicion de las tensiones en el paper

                if psiRep>0 % Material - Segun los codigos de novotny las zonas fuera del dominio no se incluyen
                    gammaSE = 1; %Valor de gamma
                


                    % Parametros de Lame
                    if ES.PEL==1
                       mu = E/(2*(1+nu));
                       lambda = nu*E/(1-nu^2);
                    elseif ES.PEL==2
                       mu = E/(2*(1+nu));
                       lambda = nu*E/( ( 1+nu)*(1-2*nu));
                    end            

                    % Tensor constitutivo
                    C = [lambda+2*mu, lambda,      0;
                        lambda,       lambda+2*mu, 0;
                        0,            0,           mu];


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

                        % Tension DE von Mises al cuadrado.
                        SVM2 = ( Sigma(1) + Sigma(2) )^2 - 3*(Sigma(1)*Sigma(2) - Sigma(3)^2);

                        % Tension relativa a la limite
                        SRelat = SVM2 / ES.SigmaM^2;

                         % Valor del integrando 
                        % para este elemento.
                        if SRelat<1e7 % Chequeo de que si es muy grande puede dar algun error esa cuenta
                            IntegElem =  gammaSE * ( ( 1 + SRelat^32 ) ^ (1/32) - 1 );
                        else
                            IntegElem = gammaSE * ( SRelat - 1 );
                        end

                        FuncTension = FuncTension + Jse * PesoG(puntoG) * Jo * IntegElem; % Se suma este punto de Gauss

                    end %En for en los puntos de Gauss
                end
            end %Fin del for de recorrer los 4 subelementos.

         % --------------------------------------------
        % Elemento material ? Segun novotny los vacios no se incluyen -----------------------------
        elseif ~ES.VA(ele)


            % Asignacion de propiedades del material estructural
            E = ES.PropMat(2,2); %Young
            nu = ES.PropMat(2,3); %Poisson


            % Obtencion de los parametros de Lame
            if ES.PEL==1
               mu = E/(2*(1+nu));
               lambda = nu*E/(1-nu^2);
            elseif ES.PEL==2
               mu = E/(2*(1+nu));
               lambda = nu*E/( ( 1+nu)*(1-2*nu));
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
             % Las tensiones se arman considerando este resultado.

            % Tensor constitutivo
            C = [lambda+2*mu, lambda,      0;
                lambda,       lambda+2*mu, 0;
                0,            0,           mu];

            Sigma= C*Eps; %Vector tension
            % En este caso dividir entre gamma es dividir entre 1 y por tanto
            % no se hace.

            % Tension DE von Mises al cuadrado.
            SVM2 = ( Sigma(1) + Sigma(2) )^2 - 3*(Sigma(1)*Sigma(2) - Sigma(3)^2);

            % Tension relativa a la limite
            SRelat = SVM2 / ES.SigmaM^2;

            % Valor del integrando 
            % para este elemento.
            % El primer 1 es para recordar que se multiplica por gamma
            % Que en el elemento material vale 1.
            if SRelat<1e7 % Chequeo de que si es muy grande puede dar algun error esa cuenta
                IntegElem = 1 * ( ( 1 + SRelat^32 ) ^ (1/32) - 1 );
            else
                IntegElem = 1 * ( SRelat - 1 );
            end


            % El chequeo de si se integra o no es en los nodos por lo que se
            % hace al final.

            % 1 PUNTO DE GAUSS EN TRIANGULOS.
            % (1 punto alcanza para integrar exactamente polinomios de
            % primer grado como es aproximadamente el caso)

            PesoG = 0.5; % Peso del punto.

            FuncTension = FuncTension + PesoG * Jo * IntegElem; % Se suma este punto de Gauss

        end % End if en la seleccion de tipo de elemento
    end % en del if de chequeo sigma
end % Fin de recorrer todos los elementos


