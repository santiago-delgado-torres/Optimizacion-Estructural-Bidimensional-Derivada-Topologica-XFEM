function ES = DerTop_Tensiones_XFEM(ES)
% Rutina para hallar la derivada topologica del funcional de tensiones
% segun el paper de Amstutz 2010 con las modificaciones hechas por Delgado
% en su tesis.
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
% ES.DTT:     Derivada topologica de las tensiones. Tiene tanta entradas como nodo la malla

% =========================================================================
% === Primeras variables ==================================================
% =========================================================================

SumVP = zeros(ES.Nnodo,1); % Se inicializa un vector de la suma de los VP.
EleCounter = zeros(ES.Nnodo,1); % Y un vector que cuenta los elementos por nodo  para calcular la media.

% =========================================================================
% === Problema Adjunto ====================================================
% =========================================================================

Uadj = ProbAdjTension(ES); % Campo de "desplazamientos" del problema adjunto.

% =========================================================================
% === Matriz de integral insoluble ========================================
% =========================================================================

nuInsol = ES.PropMat(1,3);
if ES.PEL==2
    nuInsol = nuInsol/(1-nuInsol);
end
Posiblesnu = 0:0.01:1;
indexnu = min(floor( 1+nuInsol/0.01),100); %Tengo que corregir este
nu1 = Posiblesnu(indexnu);
nu2 = Posiblesnu(indexnu+1);
InteMat1_mat2vac = load(['IntegralInsol_nu_',num2str(nu1),'_ga_1e-3.mat']); InteMat1_mat2vac = InteMat1_mat2vac.InteMat;
InteMat2_mat2vac = load(['IntegralInsol_nu_',num2str(nu2),'_ga_1e-3.mat']);InteMat2_mat2vac = InteMat2_mat2vac.InteMat;
InteMat1_vac2mat = load(['IntegralInsol_nu_',num2str(nu1),'_ga_1e3.mat']); InteMat1_vac2mat = InteMat1_vac2mat.InteMat;
InteMat2_vac2mat = load(['IntegralInsol_nu_',num2str(nu2),'_ga_1e3.mat']); InteMat2_vac2mat = InteMat2_vac2mat.InteMat;
peso = min((nuInsol-nu1)/0.01,1);
InteMat_mat2vac = InteMat1_mat2vac * (1-peso) + InteMat2_mat2vac*peso;
InteMat_vac2mat = InteMat1_vac2mat * (1-peso) + InteMat2_vac2mat*peso;

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
    
    % Funcion caracteristica del elemento
    CharFun = ES.ChequeoSigma(ele); % Esta es solo asociada a DomTen, la parte de si es estructura o no

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
        
        Uadjele = zeros(6,1); % Vector de desplazamientos nodales
        Uadjele(1:2:5) = Uadj(2*ne-1); % Se asignan los valores segun X
        Uadjele(2:2:6) = Uadj(2*ne); % Se asignan los valores segun Y        

        Uxadjele = zeros(12,1); % Vector de variables extendidas nodales
        ngx = full(ES.NodAriX(ne,ArX)); % Cuales son los grados de libertad segun x de los extendidos?
        Uxadjele(1:2:11) = Uadj(ngx); % Se asignan los valores segun x
        Uxadjele(2:2:12) = Uadj(ngx+1); % Y los segun y

        UGLAdjele = [Uadjele;Uxadjele]; % Vector de todos los GL del elemento
        
        
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
        % tercer grado como es aproximadamente el caso)

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
            
            % Ya aprovechemos y sumemos al contador:
            EleCounter(ne) = EleCounter(ne) + Jse.*ones(3,1);
            
            
            
            % Seleccion del material
            if psiRep<0 % Vacio
                E = ES.PropMat(1,2); %Young
                nu = ES.PropMat(1,3); % Poisson
                signo = -1;
                CharFunEst = 0; % Funcon caracteristica de si es elemento estructural
                gamma = 1/ES.gamma;  %El contraste es el inveso si es vacio
                b=ES.CB.Neu.Vol(1,:);  % Fuerza de volumen
                InteMat = InteMat_vac2mat;
            else
                E = ES.PropMat(2,2); %Young
                nu = ES.PropMat(2,3); % Poisson
                signo = 1;
                CharFunEst = 1; % Funcon caracteristica de si es elemento estructural
                gamma = ES.gamma;
                b=ES.CB.Neu.Vol(2,:);  % Fuerza de volumen
                InteMat = InteMat_mat2vac;
            end
            
            % Parametros de Lame
            if ES.PEL==1
               mu = E/(2*(1+nu));
               lambda = nu*E/(1-nu^2);
            elseif ES.PEL==2
               mu = E/(2*(1+nu));
               lambda = nu*E/( ( 1+nu)*(1-2*nu));
            end

            % Alpha 1 y Alpha2

            alpha2 = (lambda+3*mu) / (lambda + mu); 
            alpha1 = ( lambda + mu ) / mu; 
            

            % Tensor constitutivo
            C = [lambda+2*mu, lambda,      0;
                lambda,       lambda+2*mu, 0;
                0,            0,           mu];
            
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
                
                % Vector de "deformaciones" adjuntas elementales
                Ea = B*UGLAdjele;
                
                desplaAdj = zeros(2,1); 
                desplaAdj(1) = [ N , w1 , w2 ] * UGLAdjele(1:2:17) ; % Componente x del desplazamiento adjunto interpolado al punto de Gauss
                desplaAdj(2) =  [ N , w1 , w2 ] * UGLAdjele(2:2:18) ; % Idem componetne Y

                % Se obtiene la derivada topologica para este estado.
                DerTopG = AuxDTT(Eps,Sigma,Ea,lambda,mu,alpha1,alpha2,gamma,CharFun,ES.SigmaM,CharFunEst,signo,b,desplaAdj,InteMat);
               
                
                ValPerNodo = ( N .* DerTopG )'; % Valores por nodo en este punto
                
                
                SumVP(ne) = SumVP(ne)+ 6 * Jse * PesoG(puntoG) * ValPerNodo; % Se suma este punto de Gauss



            end %En for en los puntos de Gauss           
            

        end %Fin del for de recorrer los 4 subelementos.


    % --------------------------------------------
    % Elemento vacio -----------------------------
    elseif ES.VA(ele)
     
        CharFunEst = 0; 
        signo = -1;
        b=ES.CB.Neu.Vol(1,:);  % Fuerza de volumen
        
        % Ya aprovechemos y sumemos al contador:
        EleCounter(ne) = EleCounter(ne) + ones(3,1);
        
        
        % Asignacion de propiedades del material vacio
        E = ES.PropMat(1,2); %Young
        nu = ES.PropMat(1,3); %Poisson
        gamma = 1/ES.gamma; %El gamma del material actual (en este caso el vacio)
        
        % Parametros de Lame
        if ES.PEL==1
           mu = E/(2*(1+nu));
           lambda = nu*E/(1-nu^2);
        elseif ES.PEL==2
           mu = E/(2*(1+nu));
           lambda = nu*E/( ( 1+nu)*(1-2*nu));
        end
        
        % Alpha 1 y Alpha2
        alpha2 = (lambda+3*mu) / (lambda + mu); 
        alpha1 = ( lambda + mu ) / mu;
        
        
        % Como en el elemento tradicional tanto las tensiones, como las
        % deformaciones e incluso el Ea del problema adjunto son UNIFORMES.
        % Alcanza calcularlas en este elemento sin tener que hacer nada de
        % puntos de Gauss o similar.
        
        Uele = zeros(6,1); % Vector de desplazamientos nodales
        Uele(1:2:5) = ES.U(2*ne-1); % Se asignan los valores segun X
        Uele(2:2:6) = ES.U(2*ne); % Se asignan los valores segun Y
        
        Uadjele = zeros(6,1); %Vector de "desplazamientos" adjuntos nodales
        Uadjele(1:2:5) = Uadj(2*ne-1); % Se asignan los valores segun X
        Uadjele(2:2:6) = Uadj(2*ne); % Se asignan los valores segun y
        

        % Matriz relacion desplazamiento-deformacion
        B = zeros(3,6);
        B(1,1:2:5)=dN_dxy(1,:);
        B(2,2:2:6)=dN_dxy(2,:);
        B(3,2:2:6)=dN_dxy(1,:);
        B(3,1:2:5)=dN_dxy(2,:);

        % Vector de deformaciones elementales
        Eps = B*Uele;
        
        % Vector de "deformaciones" adjuntas elementales.
        Ea = B*Uadjele;

        % Tensor constitutivo
        C = [lambda+2*mu, lambda,      0;
            lambda,       lambda+2*mu, 0;
            0,            0,           mu];

        % Vector de tensiones elementales
        Sigma=C*Eps;
        
        % 1 PUNTO DE GAUSS EN TRIANGULOS.
        % (1 punto alcanza para integrar exactamente polinomios de
        % primer grado como es aproximadamente el caso)

        N1 = 1/3; % Funcion de forma N1 clasica para el elemento en el punto de Gauss
        N2 = 1/3; % Idem N2
        N3 = 1/3; % Idem N3
        PesoG = 0.5; % Peso del punto.
        N=[N1,N2,N3]; % En formato vector estas 3
        
        desplaAdj = zeros(2,1);
        desplaAdj(1) = N*Uadjele(1:2:5);
        desplaAdj(2) = N*Uadjele(2:2:6);
        
        % Se obtiene la derivada topologica para este estado.
        DerTopG = AuxDTT(Eps,Sigma,Ea,lambda,mu,alpha1,alpha2,gamma,CharFun,ES.SigmaM,CharFunEst,signo,b,desplaAdj,InteMat_vac2mat);
               
        
        
     
        ValPerNodo = ( N .* DerTopG)';
       
        SumVP(ne) = SumVP(ne)+ 6 * PesoG * ValPerNodo; % Se suma este punto de Gauss
     
     % --------------------------------------------
     % Elemento tradicional -----------------------
    else
             
        % Ya aprovechemos y sumemos al contador:
        EleCounter(ne) = EleCounter(ne) + ones(3,1);
        
        
        % Asignacion de propiedades del material estructural
        E = ES.PropMat(2,2); %Young
        nu = ES.PropMat(2,3); %Poisson
        gamma = ES.gamma;
        signo = 1;
        CharFunEst = 1;
        b=ES.CB.Neu.Vol(2,:);  % Fuerza de volumen
        
        % Parametros de Lame
        if ES.PEL==1
           mu = E/(2*(1+nu));
           lambda = nu*E/(1-nu^2);
        elseif ES.PEL==2
           mu = E/(2*(1+nu));
           lambda = nu*E/( ( 1+nu)*(1-2*nu));
        end
        % Alpha 1 y Alpha2
        alpha2 = (lambda+3*mu) / (lambda + mu); 
        alpha1 = ( lambda + mu ) / mu;
        
        % Como en el elemento tradicional tanto las tensiones, como las
        % deformaciones e incluso el Ea del problema adjunto son UNIFORMES.
        % Alcanza calcularlas en este elemento sin tener que hacer nada de
        % puntos de Gauss o similar.
        
        Uele = zeros(6,1); % Vector de desplazamientos nodales
        Uele(1:2:5) = ES.U(2*ne-1); % Se asignan los valores segun X
        Uele(2:2:6) = ES.U(2*ne); % Se asignan los valores segun Y
        
        Uadjele = zeros(6,1); %Vector de "desplazamientos" adjuntos nodales
        Uadjele(1:2:5) = Uadj(2*ne-1); % Se asignan los valores segun X
        Uadjele(2:2:6) = Uadj(2*ne); % Se asignan los valores segun y
        

        % Matriz relacion desplazamiento-deformacion
        B = zeros(3,6);
        B(1,1:2:5)=dN_dxy(1,:);
        B(2,2:2:6)=dN_dxy(2,:);
        B(3,2:2:6)=dN_dxy(1,:);
        B(3,1:2:5)=dN_dxy(2,:);

        % Vector de deformaciones elementales
        Eps = B*Uele;
        
        % Vector de "deformaciones" adjuntas elementales.
        Ea = B*Uadjele;

        % Tensor constitutivo
        C = [lambda+2*mu, lambda,      0;
            lambda,       lambda+2*mu, 0;
            0,            0,           mu];

        % Vector de tensiones elementales
        Sigma=C*Eps; 
        
        
        % 1 PUNTO DE GAUSS EN TRIANGULOS.
        % (1 punto alcanza para integrar exactamente polinomios de
        % primer grado como es aproximadamente el caso)

        N1 = 1/3; % Funcion de forma N1 clasica para el elemento en el punto de Gauss
        N2 = 1/3; % Idem N2
        N3 = 1/3; % Idem N3
        PesoG = 0.5; % Peso del punto.
        N=[N1,N2,N3]; % En formato vector estas 3
        
        desplaAdj = zeros(2,1);
        desplaAdj(1) = N*Uadjele(1:2:5);
        desplaAdj(2) = N*Uadjele(2:2:6);
        
        
        % Se obtiene la derivada topologica para este estado.
        DerTopG = AuxDTT(Eps,Sigma,Ea,lambda,mu,alpha1,alpha2,gamma,CharFun,ES.SigmaM,CharFunEst,signo,b,desplaAdj,InteMat_mat2vac);
               
        
     
        ValPerNodo = ( N .* DerTopG)';
       
        SumVP(ne) = SumVP(ne)+ 6 * PesoG * ValPerNodo; % Se suma este punto de Gauss
        

    end % End if en la seleccion de tipo de elemento

end % Fin de recorrer todos los elementos


ES.DTT = SumVP./EleCounter; % Se calcula la media
