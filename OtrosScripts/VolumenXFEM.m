function Omega = VolumenXFEM(ES)
  
% Funcion que calcula el volumen de material estructural en la malla
%
% Devuelve Omega: Volumen actual;

% =========================================================================
% === Iteracion en los elementos ==========================================
% =========================================================================

Omega=0; % Comienza como nulo


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
        
         % FUNCIONES DE INTERPOLACION CLASICAS
        N1 = @(eta,xi) 1-eta-xi;
        N2 = @(eta,xi) eta;
        N3 = @(eta,xi) xi;
    
        N = @(eta,xi) [N1(eta,xi),N2(eta,xi),N3(eta,xi)]; % Vector de las funciones de forma en formato funcion.
        
        psie = ES.psi(ne); %Funcion de nivel para los nodos del elemento.
        
        % PUNTOS INTERMEDIOS PARA DIVISION EN SUBELEMENTOS PARA INTEGRACION
        
        % Punto entre nodo 1 y 2. En este xi = 0
        if psie(1)*psie(2)<0 
            etaR1 = psie(1)/(psie(1)-psie(2));
        else
            etaR1 = 0.5;
        end
        
        % Punto entre nodo 2 y 3. En este xi = 1-eta
        if psie(2)*psie(3)<0
            etaR2 = psie(3)/(psie(3)-psie(2));
        else
            etaR2 = 0.5;
        end  
        
        % Punto entre nodo 3 y 1. En este eta=0.        
        if psie(3)*psie(1)<0
            xiR3 = psie(1)/(psie(1)-psie(3));
        else
            xiR3 = 0.5;
        end  
        
        % 1 PUNTOS DE GAUSS EN TRIANGULOS.
        % (1 punto alcanza para integrar exactamente polinomios de
        % ctes como es el caso)
        % Se integra 1 si psiRep>0 y 0 si no
        
        PesoG = 0.5; % Peso de los puntos.
        
       
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
        
        for subEl = 1:4 % Se recorren los 4 subelementos
            
            etase = etaSubElements(:,subEl); % Coordenadas eta de este subtriangulo.
            xise = xiSubElements(:,subEl); % Idem para la coordenada xi.
            
            Jse = det( dN_detachi*[etase xise] ); % Determinante jacobiano de este subtriangulo en el triangulo intrinseco
            
            % Siguiente for es para obtener los valores de psi en los
            % vertices del subtriangulo.
            psise = zeros(3,1); % Se inicializa como nulo
            for punto = 1:3
                psise(punto) = N( etase(punto) , xise(punto) ) * psie; %Se obtiene como interpolacion lineal del otro triangulo.
            end
            
            psiRep=mean(psise); % Por la subdivision los psi de los vertices de los subtriangulos, tienen todos el mismo signo o son nulos.
            % Por ende un valor representativo para ver el material de
            % este subtriangulo en particular es calcular la media.
            
            % Suma de volumen
            
            if psiRep>0
              Omega = Omega + PesoG * Jo * Jse;
            end % Si psiRep<0 no se hace nada
            
            
            
        end %Fin del for de recorrer los 4 subelementos.
         
       
    % --------------------------------------------
    % Elemento material -----------------------------
    elseif ~ES.VA(ele)
    
        PesoG=0.5;
        
        Omega = Omega + Jo * PesoG; 
    
        
    end % End if en la seleccion de tipo de elemento. 
    % No se incluye un caso de elemento vacio pues esos son de volumen nulo
end %Fin de la iteracion en los elementos.
