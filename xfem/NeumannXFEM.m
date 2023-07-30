function ES = NeumannXFEM(ES)

% Rutina para hallar el vector de fuerzas externas previo a los efectos de
% desplazamientos impuestos.
%
% Devuelve o actualiza:
%
% Fext: Vector de fuerzas externas

% =========================================================================
% === Inicializacion del vector ===========================================
% =========================================================================

ES.Fext = zeros(2*ES.Nnodo + ES.NGLX,1); 

% =========================================================================
% === Fuerzas nodales =====================================================
% =========================================================================

for i = 1:size(ES.CB.Neu.Punt,1) % Recorremos todas las entradas de la matriz asociada
    
    ES.Fext( 2*ES.CB.Neu.Punt(i,1) - 1 ) = ES.Fext( 2*ES.CB.Neu.Punt(i,1) -1 ) + ES.CB.Neu.Punt(i,2) ; % Se suma la componente segun X
    
    ES.Fext( 2*ES.CB.Neu.Punt(i,1)  ) = ES.Fext( 2*ES.CB.Neu.Punt(i,1)  ) + ES.CB.Neu.Punt(i,3) ; % Se suma la componente segun Y
    
end

% =========================================================================
% === Fuerzas en Bordes ===================================================
% =========================================================================

for i = 1:size(ES.CB.Neu.Bord,1) % Se recorren todas las entradas de la matriz asociada.
    
    % Pequeño chequeo de que la definicion del tipo de fuerza sea correcta

    if ( ES.CB.Neu.Bord(i,1) ~= 1 ) && ( ES.CB.Neu.Bord(i,1) ~= 2 )
        error(['ES.CB.Neu.Bord desconoce el tipo ',num2str(ES.CB.Neu.Bord(i,1)) ,' que se dio en la entrada ',num2str(i),'. La misma debe valer 1 (Segun X) o 2 (Segun Y)'])
    end

    % Coeficiente para asignar si la fuerza es segun x o segun y.
    if ES.CB.Neu.Bord(i,1) == 1 
        Coeficiente=-1;
    else
        Coeficiente=0;
    end
    % la logica es luego en la asignacion es 2*numero de nodo + Coeficiente
    % Para los GL no extendidos.
    % En los GL extendidos es NodAriX + Coeficiente + 1
    
    ele = ES.CB.Neu.Bord(i,2) ; % Variable auxiliar con el numero de elemento para no tener que volver a escribir todo
    
    ne  = ES.Melem(ele,[3:5,3]); % Los nodos del elemento. Lo de repetir el 3 es para un truquito con el siguiente paso
    
    nb = ne( ES.CB.Neu.Bord(i,3):(ES.CB.Neu.Bord(i,3)+1) ); %Nodos del borde donde esta la fuerza
    
    Jo = sqrt ( ( ES.Mnodo(nb(1),2) - ES.Mnodo(nb(2),2) ).^2 + ( ES.Mnodo(nb(1),3) - ES.Mnodo(nb(2),3) ).^2 ) / 2; % Determinante jacobiano para integrar 
    
    % Incorporacion a los grados de libertad no extendidos
    % Al ser una fuerza constante, y la funcion de forma ser lineal,
    % alcanza con un unico punto de Gauss (De peso 2 y ubicado en eta=0)
    
    ES.Fext( 2*nb + Coeficiente ) = ES.Fext( 2*nb + Coeficiente ) + Jo * ES.CB.Neu.Bord(i,4); 
    % Faltaria multiplicar el agregado por el peso de Gauss (2) y por el
    % valor de la funcion de forma en el punto (0.5). Como este producto es
    % 1 se obvia.
    
    % Incorporacion a los grados de libertad extendidos si corresponde.
    
    if ES.psi(nb(1))*ES.psi(nb(2))<0 % Es decir, si ademas de extendido, se da que en este borde en particular halla interfase
        % Esto importa porque si en el borde de interes no hay interfase,
        % ambas funciones "g" son nulas en el borde y por ende tambien las
        % funciones de forma wi. 
        % Siendo estas nulas, el trabajo virtual de la fuerza externa en el
        % borde es nula para cualquier valor del S virtual y por ende la
        % componente en el vector de fuerzas lo es.
        
        ArX = ES.EGLX(ele,:); %Numero de aristas extendidas de este elemento
        NArX = ES.AriX(ArX,:); %Numero de nodo de esas aristas
        if ( nb(1)==NArX(1,1) || nb(1)==NArX(1,2) ) && ( nb(2)==NArX(1,1) || nb(2) == NArX(1,2) ) % La fuerza es en la arista extendida 1 
            AriX = 1; 
        else % Es la arista 2
            AriX=2;
        end
        
            
        % Definiendo como eta = -1 el punto en nb(1) y eta = 1 el punto en
        % nb(2). El eta donde esta la interfas es:
        
        etaint = ( 2*ES.psi(nb(1)) / ( ES.psi(nb(1)) - ES.psi(nb(2)) ) ) -1;
        
        % La logica es integrar en dos subdominios, de forma analoga a lo
        % hecho para la matriz de rigidez.
        
        % El valor de g en la interfaz es:
        ginte = 1; %Solo nos interesa considerar la que vale 1 pues la otra arista es toda nula en el borde
        
        GL1 = ES.NodAriX(nb(1),ArX(AriX)) + Coeficiente + 1; % Grado de libertad asociado al primer nodo
        GL2 = ES.NodAriX(nb(2),ArX(AriX)) + Coeficiente + 1; %Idem para el segundo nodo
        % El tercer nodo no le importa pues N3 vale 0 en todo el borde
        
        Jse = [( etaint + 1 ) / 2 ; % Determinante jacobiano del primer tramo
              ( 1 - etaint ) / 2 ]; % Determinante jacobiano del segundo tramo
               
        % Se integra mediante 2 puntos de gauss ya que la funcion de interpolacion es cuadratica
        % Estos puntos son en \pm 1/sqrt(3) y su peso es 1
        pesoG = 1;
        
        etaPG = [-1/sqrt(3) 1/sqrt(3) ]; 
        
        for subtr = 1:2 % Se recorren los 2 subtramos 
            for PG = 1:2 %Y los 2 puntos de Gauss
        
                if subtr==1 
                    gPG =  ginte * (etaPG(PG)+1) / 2 ;%Valor de g en el punto de integracion de Gauss si es el primer subtramo
                    etaeq = -1 + (etaint+1) * ( etaPG(PG) + 1 )/2; %Valor de eta (del borde completo) en el punto de Gauss.
                else 
                    gPG = ginte * ( 1 - etaPG(PG) ) /2 ; %Valor de g en el punto de integracion de Gauss si es el segundo subtramo
                    etaeq = 1 + (etaint-1) * ( 1-etaPG(PG)) / 2; %Valor de eta del borde completo en el punto de Gauss.
                end
                               

                N1PG = (1-etaeq) / 2 ; %Funcion de forma para el primer nodo en el Punto de integracion de Gauss.

                N2PG = (etaeq+1) / 2; % Idem pero para el segundo nodo
        
                w1PG = gPG * N1PG ; % Funcion de interpolacion para el grado de libertad extendido del primer nodo en el borde para este punto
                w2PG = gPG * N2PG; % Idem para el segundo nodo.
                
                ES.Fext( GL1 ) = ES.Fext( GL1 ) + Jo * Jse(subtr) * pesoG * w1PG  * ES.CB.Neu.Bord(i,4); % Se agrega para el primer GL
                ES.Fext( GL2 ) = ES.Fext( GL2 ) + Jo * Jse(subtr) * pesoG * w2PG  * ES.CB.Neu.Bord(i,4); % Se agrega para el segundo GL
            end %En recorrer puntos de gauss
        end % En recorrer subtramos     
    end %Del if de que pasa si hay que incluir en los extendidos
end % En recorrer cargas en bordes

% =========================================================================
% === Fuerzas de Peso Propio ==============================================
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
    
    % --------------------------------------------
    % Elemento extendido -------------------------
    if ES.EI(ele)
        % De manera analoga a la matriz de rigidez tangente hay que subdividir
        % en subtriangulos para hacer la integral. 
        % Aun en el caso de las componentes tradicionales ya que el b 
        % es potencialmente distinto en la zona de material vacio y la zona de material
        % estructural.
        
                
        %Funcion de nivel para los nodos del elemento.
        psie = ES.psi(ne); 
        
         
        ArX = ES.EGLX(ele,:); %Numero de aristas extendidas de este elemento
        NArX = ES.AriX(ArX,:); %Numero de nodo de esas aristas
                
        % PUNTOS INTERMEDIOS PARA DIVISION EN SUBELEMENTOS PARA INTEGRACION
        
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
            
            Jse = det( dN_detachi*[etase xise] ); % Determinante jacobiano de este subtriangulo en el triangulo intrinseco
                                
             % Obtencion del material
            % Se obtiene cuanto vale el level set en un punto intermedio
            % del subelemento
            etaeval =  [1/3 1/3 1/3]* etase; % Interpolacion para el valor de eta.
            xieval = [1/3 1/3 1/3] *xise; %Idem para el valor de xi
            
            psiRep = [1-etaeval-xieval,etaeval,xieval]*psie; % Valor del level set representativo del SE
            
            % Seleccion de la fuerza de volumen
            
            if psiRep<0 % MAterial vacio
                b=ES.CB.Neu.Vol(1,:); 
            else % Material estructural
                b=ES.CB.Neu.Vol(2,:); 
            end
            
            % Integracion en puntos de Gauss
            
            for puntoG=1:3 % Para los 3 puntos de Gauss anteriormente definidos
                
                % Los puntos de Gauss son en sistemas de coordenadas
                % intrinsecos a los subtriangulos.
                
                % Para saber el valor de N/w se debe hayar los valores de eta
                % y xi en el triangulo completo:
                N1se = 1-etaG(puntoG)-xiG(puntoG) ; % Funcion de interpolacion clasica N1 en el subtriangulo
                N2se = etaG(puntoG); %Idem funcion N2
                N3se = xiG(puntoG); % Idem funcion N3
                Nse = [N1se,N2se,N3se]; % Vector de estas 3 funciones
                etaeval = Nse * etase; % Interpolacion para el valor de eta.
                xieval = Nse *xise; %Idem para el valor de xi
                 
                % Funciones de interpolacion clasicas en este punto:
                N1 = 1-etaeval-xieval; % funcion de interpolacion clasica N1 para el elemento completo en este punto
                N2 = etaeval; %Idem N2
                N3 = xieval; %Idem N3
                N = [N1,N2,N3]; % vector de las mismas.    
    
                % Funciones para definir los w
                g1 = Nse * Ar1;
                g2 = Nse * Ar2;
                w1 = N*g1;
                w2 = N*g2;
                
                % Fuerzas para este punto de Gauss
                
                Ftx = Jo * Jse * PesoG(puntoG) * N * ES.esp * b(1); % Componente segun x en los GL tradicionales
                Fty = Jo * Jse * PesoG(puntoG) * N * ES.esp * b(2); % Componente segun y en los GL tradicionales
                Fxx1 = Jo * Jse * PesoG(puntoG) * w1 * ES.esp * b(1); % Componente segun x en los GL extendidos
                Fxy1 = Jo * Jse * PesoG(puntoG) * w1 * ES.esp * b(2); % Componente segun y en los GL extendidos
                Fxx2 = Jo * Jse * PesoG(puntoG) * w2 * ES.esp * b(1); % Componente segun x en los GL extendidos
                Fxy2 = Jo * Jse * PesoG(puntoG) * w2 * ES.esp * b(2); % Componente segun y en los GL extendidos
               
                % Incorporacion al vector de fuerzas
       
                ES.Fext( 2*ne - 1 ) = ES.Fext( 2*ne - 1 ) + Ftx'; % Agregado a los Grados de libertad tradicionales s/ x
                ES.Fext( 2*ne     ) = ES.Fext( 2*ne     ) + Fty'; % Agregado a los GL tradicionales s/ y
                ES.Fext( full( ES.NodAriX(ne,ArX(1)) ) ) = ES.Fext( full( ES.NodAriX(ne,ArX(1)) ) ) + Fxx1'; % Agregado a los GL extendidos segun x
                ES.Fext( full( ES.NodAriX(ne,ArX(1)) )+1 ) = ES.Fext( full( ES.NodAriX(ne,ArX(1)) )+1 ) + Fxy1'; % Agregado a los GL extendidos segun y
                ES.Fext( full( ES.NodAriX(ne,ArX(2)) ) ) = ES.Fext( full( ES.NodAriX(ne,ArX(2)) ) ) + Fxx2'; % Agregado a los GL extendidos segun x
                ES.Fext( full( ES.NodAriX(ne,ArX(2)) )+1 ) = ES.Fext( full( ES.NodAriX(ne,ArX(2)) )+1 ) + Fxy2'; % Agregado a los GL extendidos segun y

            end %En for en los puntos de Gauss
            
        end %Fin del for de recorrer los 4 subelementos.        
        
        
    
    % --------------------------------------------
    % Elemento vacio -----------------------------
    elseif ES.VA(ele)
        b=ES.CB.Neu.Vol(1,:); %Valores de la fuerza de volumen segun x y segun y para el material vacio
        
        % Al ser las funciones de forma lineales, alcanza con integrar en
        % un unico punto de Gauss. 
        % Este punto tiene peso 0,5 y los valores son etaG=xiG = 1/3;
        % En este caso se sabe que las 3 funciones de forma valen 1/3.
        
        PesoG = 0.5; % Peso de Gauss
        N = 1/3; %Valor de las funciones de forma en el punto de Gauss
        
        Fx = Jo * PesoG * N * ES.esp * b(1); %Fuerza a agregar segun x
        Fy = Jo * PesoG * N * ES.esp * b(2);  % Fuerza a agregar segun y
        
        ES.Fext( 2*ne - 1 ) = ES.Fext( 2*ne - 1 ) + Fx*ones(3,1); %Agregados a los grados de libertad correspondientes
        ES.Fext( 2*ne     ) = ES.Fext( 2*ne     ) + Fy*ones(3,1); %Idem
    
    % --------------------------------------------
    % Elemento tradicional -----------------------    
    else
        
        b=ES.CB.Neu.Vol(2,:); %Valores de la fuerza de volumen segun x y segun y para el material tradicional
        
        % Al ser las funciones de forma lineales, alcanza con integrar en
        % un unico punto de Gauss. 
        % Este punto tiene peso 0,5 y los valores son etaG=xiG = 1/3;
        % En este caso se sabe que las 3 funciones de forma valen 1/3.
        
        PesoG = 0.5; % Peso de Gauss
        N = 1/3; %Valor de las funciones de forma en el punto de Gauss
        
        Fx = Jo * PesoG * N * ES.esp * b(1); %Fuerza a agregar segun x
        Fy = Jo * PesoG * N * ES.esp * b(2);  % Fuerza a agregar segun y
        
        ES.Fext( 2*ne - 1 ) = ES.Fext( 2*ne - 1 ) + Fx*ones(3,1); %Agregados a los grados de libertad correspondientes
        ES.Fext( 2*ne     ) = ES.Fext( 2*ne     ) + Fy*ones(3,1); %Idem
        
        
    end % Del if de que tipo de elemento es.
    

end % En recorrer todos los elementos







end 