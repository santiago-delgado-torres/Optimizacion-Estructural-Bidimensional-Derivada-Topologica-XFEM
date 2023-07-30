function [] = figXFEM(ES, tipo, vacum, bord, Scale)

% Rutina para graficar la malla de elementos finitos

% -------------------------------------------------------------------------
% Datos -------------------------------------------------------------------
% Vacum es un logico, si es 0 no se muestran los elementos (o subtriangulos)
% vacios
% Bord es otro logico. Si es 0 no se muestran los bordes de la malla
% Scale es la escala para la deformada unicamente


if strcmp(tipo,'UX') %Desplazamientos segun X
  
    Mnodo = [ES.Mnodo(:,2:3); zeros(sum(ES.EI)*3,2)]; % Se inicia matriz de nodos agregando 3 filas por cada elemento extendido
    Melem = [ES.Melem(:,3:5); zeros(sum(ES.EI)*3,3)]; % Idem para la matriz de elementos
    % Como la idea es agregar algunos mas a medida que se encuentren interfases
    % Estos son indices actuales de la posicion ingresada de nodo
    % y de elemento respectivamente
    PosN = ES.Nnodo; 
    PosE = ES.Nelem;
    Ui = [ES.U(1:2:(2*ES.Nnodo-1)); zeros(sum(ES.EI)*3,1)]; % Se inicia vector de la variable a pintar
    
    Vac = zeros(size(Melem,1),1); % Inicializado vector logico para chequeo de vacios
    
    for ele = 1:ES.Nelem
        % Nodos
        ne = ES.Melem(ele,3:5);
        
        % Geometria del elemento
        Xe = ES.Mnodo(ne,2);
        Ye = ES.Mnodo(ne,3);
        XYe=[Xe,Ye];
        
        % Funcion de nivel en el elemento
        psie = ES.psi(ne);
        
        
        if ~ES.EI(ele) && mean(psie)<0 % Si no hay interfaz es facil chequear si es vacio o material
          % Ademas, si era material sin interfaz no hay que hacer nada en Vac,
            Vac(ele) = 1;
          
        end
        
        % ---------------------------------------------------------------------
        % Elemento con interfase ----------------------------------------------
        if ES.EI(ele)
            
            ArX = ES.EGLX(ele,:); %Numero de aristas extendidas de este elemento
            NArX = ES.AriX(ArX,:); %Numero de nodo de esas aristas

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
            xiR1=0;

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
            xiR2=1-etaR2;

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
            etaR3=0;
            
            % Nodos adicionales
            N = zeros(3,3);
            N(1,1) = 1-etaR1-xiR1; N(1,2) = etaR1; N(1,3) = xiR1;
            N(2,1) = 1-etaR2-xiR2; N(2,2) = etaR2; N(2,3) = xiR2;
            N(3,1) = 1-etaR3-xiR3; N(3,2) = etaR3; N(3,3) = xiR3;
            XYI = N*XYe; % Coordenadas de los 3 nodos adicionales
            Mnodo(PosN+1:PosN+3,:) = XYI;
            Melem(ele,:) = [ne(1), PosN+1, PosN+3]; % sustituyo elemento
            % Este es entre el primer nodo, el punto intermedio entre 1 y 2 y el intermedio entre 1 y3
            
            Melem(PosE+1,:) = [PosN+1, ne(2), PosN+2]; % Este en el punto entre 1 y 2, el nodo 2 y el intermedio entre 2 y 3
            Melem(PosE+2,:) = [PosN+1, PosN+2, PosN+3]; % Este es el central
            Melem(PosE+3,:) = [PosN+3, PosN+2, ne(3)]; %El que queda

            
            % Chequeo de vacios
            PsiABC = N*psie; % Valores de psi para los 3 nodos adicionales
            
            % Elemento sustituido
            if (psie(1) + PsiABC(1) + PsiABC(3))<0
              Vac(ele) = 1;
            end
            % Primer elemento nuevo
            if (PsiABC(1) + psie(2) + PsiABC(2))<0
              Vac(PosE+1)=1;
            end
            % Segundo
            if (PsiABC(1) + PsiABC(2)+PsiABC(3))<0
              Vac(PosE+2)=1;
            end
            % Tercero
            if (PsiABC(3) + PsiABC(2)+psie(3))<0
              Vac(PosE+3)=1;
            end
            
            
            PosE = PosE+3;
            
            % Eta y Xi elemento interfase
            ETI = [etaR1; etaR2; etaR3];
            XII = [xiR1; xiR2; xiR3];
            
            % Valores de las aristas en elemnto interfase
            AR1 = [Ar1R1;Ar1R2;Ar1R3];
            AR2 = [Ar2R1;Ar2R2;Ar2R3];
            
            N = zeros(1,9);
            
            % Valores nodales
            UX = [ES.U(2*ne-1); ES.U( reshape( ES.NodAriX(ne,ArX) , [6,1]) )]; 

            for nod = 1:3
                % Coeficientes de la funcion de forma de nodo extendido
                et = ETI(nod); xi = XII(nod);
                N(1) = 1-et-xi; N(2) = et; N(3) = xi;
                
                % Calculo de las g
                g1 = AR1(nod); %g vale simplemente la auxiliar arista en elemento inferfase
                g2 = AR2(nod);
                N(4:6)= N(1:3)*g1;
                N(7:9) = N(1:3)*g2;
                
                UI = N*UX;
                
                Ui(PosN+nod) = UI;
            end
            
            PosN = PosN + 3;
        end
    end
   if vacum==0 % En caso de querer borrar los elementos vacios para el plot
     Melem(Vac==1,:)=[]; % Se los elimina de la matriz de elementos para el trisurf
   end  
    % Escribire en la figura abierta
    
   if bord==0
      Opcion='none';
    else
      Opcion='k';
    end
    
    hold on;
    trisurf(Melem,Mnodo(:,1),Mnodo(:,2),Ui,Ui,'FaceColor','interp','EdgeColor',Opcion);
    
elseif strcmp(tipo,'UY') %Desplazamientos segun Y
  
    Mnodo = [ES.Mnodo(:,2:3); zeros(sum(ES.EI)*3,2)]; % Se inicia matriz de nodos agregando 3 filas por cada elemento extendido
    Melem = [ES.Melem(:,3:5); zeros(sum(ES.EI)*3,3)]; % Idem para la matriz de elementos
    % Como la idea es agregar algunos mas a medida que se encuentren interfases
    % Estos son indices actuales de la posicion ingresada de nodo
    % y de elemento respectivamente
    PosN = ES.Nnodo; 
    PosE = ES.Nelem;
    Ui = [ES.U(2:2:(2*ES.Nnodo)); zeros(sum(ES.EI)*3,1)]; % Se inicia vector de la variable a pintar
    
    Vac = zeros(size(Melem,1),1); % Inicializado vector logico para chequeo de vacios
    
    for ele = 1:ES.Nelem
        % Nodos
        ne = ES.Melem(ele,3:5);
        
        % Geometria del elemento
        Xe = ES.Mnodo(ne,2);
        Ye = ES.Mnodo(ne,3);
        XYe=[Xe,Ye];
        
        % Funcion de nivel en el elemento
        psie = ES.psi(ne);
        
        
        if ~ES.EI(ele) && mean(psie)<0 % Si no hay interfaz es facil chequear si es vacio o material
          % Ademas, si era material sin interfaz no hay que hacer nada en Vac,
            Vac(ele) = 1;
          
        end
        
        % ---------------------------------------------------------------------
        % Elemento con interfase ----------------------------------------------
        if ES.EI(ele)
            
            ArX = ES.EGLX(ele,:); %Numero de aristas extendidas de este elemento
            NArX = ES.AriX(ArX,:); %Numero de nodo de esas aristas

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
            xiR1=0;

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
            xiR2=1-etaR2;

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
            etaR3=0;
            
            % Nodos adicionales
            N = zeros(3,3);
            N(1,1) = 1-etaR1-xiR1; N(1,2) = etaR1; N(1,3) = xiR1;
            N(2,1) = 1-etaR2-xiR2; N(2,2) = etaR2; N(2,3) = xiR2;
            N(3,1) = 1-etaR3-xiR3; N(3,2) = etaR3; N(3,3) = xiR3;
            XYI = N*XYe; % Coordenadas de los 3 nodos adicionales
            Mnodo(PosN+1:PosN+3,:) = XYI;
            Melem(ele,:) = [ne(1), PosN+1, PosN+3]; % sustituyo elemento
            % Este es entre el primer nodo, el punto intermedio entre 1 y 2 y el intermedio entre 1 y3
            
            Melem(PosE+1,:) = [PosN+1, ne(2), PosN+2]; % Este en el punto entre 1 y 2, el nodo 2 y el intermedio entre 2 y 3
            Melem(PosE+2,:) = [PosN+1, PosN+2, PosN+3]; % Este es el central
            Melem(PosE+3,:) = [PosN+3, PosN+2, ne(3)]; %El que queda

            
            % Chequeo de vacios
            PsiABC = N*psie; % Valores de psi para los 3 nodos adicionales
            
            % Elemento sustituido
            if (psie(1) + PsiABC(1) + PsiABC(3))<0
              Vac(ele) = 1;
            end
            % Primer elemento nuevo
            if (PsiABC(1) + psie(2) + PsiABC(2))<0
              Vac(PosE+1)=1;
            end
            % Segundo
            if (PsiABC(1) + PsiABC(2)+PsiABC(3))<0
              Vac(PosE+2)=1;
            end
            % Tercero
            if (PsiABC(3) + PsiABC(2)+psie(3))<0
              Vac(PosE+3)=1;
            end
            
            
            PosE = PosE+3;
            
            % Eta y Xi elemento interfase
            ETI = [etaR1; etaR2; etaR3];
            XII = [xiR1; xiR2; xiR3];
            
            % Valores de las aristas en elemnto interfase
            AR1 = [Ar1R1;Ar1R2;Ar1R3];
            AR2 = [Ar2R1;Ar2R2;Ar2R3];
            
            N = zeros(1,9);
            
            % Valores nodales
            UX = [ES.U(2*ne); ES.U( reshape( ES.NodAriX(ne,ArX) +1 , [6,1]) )]; 

            for nod = 1:3
                % Coeficientes de la funcion de forma de nodo extendido
                et = ETI(nod); xi = XII(nod);
                N(1) = 1-et-xi; N(2) = et; N(3) = xi;
                
                % Calculo de las g
                g1 = AR1(nod); %g vale simplemente la auxiliar arista en elemento inferfase
                g2 = AR2(nod);
                N(4:6)= N(1:3)*g1;
                N(7:9) = N(1:3)*g2;
                
                UI = N*UX;
                
                Ui(PosN+nod) = UI;
            end
            
            PosN = PosN + 3;
        end
    end
   if vacum==0 % En caso de querer borrar los elementos vacios para el plot
     Melem(Vac==1,:)=[]; % Se los elimina de la matriz de elementos para el trisurf
   end  
    % Escribire en la figura abierta
    
   if bord==0
      Opcion='none';
    else
      Opcion='k';
    end
    
    hold on;
    trisurf(Melem,Mnodo(:,1),Mnodo(:,2),Ui,Ui,'FaceColor','interp','EdgeColor',Opcion);
    
elseif strcmp(tipo,'DEF') %Deformada
  
    Mnodo = [ES.Mnodo(:,2:3); zeros(sum(ES.EI)*3,2)]; % Se inicia matriz de nodos agregando 3 filas por cada elemento extendido
    Melem = [ES.Melem(:,3:5); zeros(sum(ES.EI)*3,3)]; % Idem para la matriz de elementos
    % Como la idea es agregar algunos mas a medida que se encuentren interfases
    % Estos son indices actuales de la posicion ingresada de nodo
    % y de elemento respectivamente
    PosN = ES.Nnodo; 
    PosE = ES.Nelem;
    UiX = [ES.U(1:2:(2*ES.Nnodo-1)); zeros(sum(ES.EI)*3,1)]; % Se inicia vector de la variable con los desplazamientos X
    UiY = [ES.U(2:2:(2*ES.Nnodo)); zeros(sum(ES.EI)*3,1)]; % Se inicia vector de la variable con los desplazamientos Y
    
    Vac = zeros(size(Melem,1),1); % Inicializado vector logico para chequeo de vacios
    
    
    for ele = 1:ES.Nelem
        % Nodos
        ne = ES.Melem(ele,3:5);
        
        % Geometria del elemento
        Xe = ES.Mnodo(ne,2);
        Ye = ES.Mnodo(ne,3);
        XYe=[Xe,Ye];
        
        % Funcion de nivel en el elemento
        psie = ES.psi(ne);
        
        
        if ~ES.EI(ele) && mean(psie)<0 % Si no hay interfaz es facil chequear si es vacio o material
          % Ademas, si era material sin interfaz no hay que hacer nada en Vac,
            Vac(ele) = 1;
          
        end
        
        % ---------------------------------------------------------------------
        % Elemento con interfase ----------------------------------------------
        if ES.EI(ele)
            
            ArX = ES.EGLX(ele,:); %Numero de aristas extendidas de este elemento
            NArX = ES.AriX(ArX,:); %Numero de nodo de esas aristas

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
            xiR1=0;

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
            xiR2=1-etaR2;

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
            etaR3=0;
            
            % Nodos adicionales
            N = zeros(3,3);
            N(1,1) = 1-etaR1-xiR1; N(1,2) = etaR1; N(1,3) = xiR1;
            N(2,1) = 1-etaR2-xiR2; N(2,2) = etaR2; N(2,3) = xiR2;
            N(3,1) = 1-etaR3-xiR3; N(3,2) = etaR3; N(3,3) = xiR3;
            XYI = N*XYe; % Coordenadas de los 3 nodos adicionales
            Mnodo(PosN+1:PosN+3,:) = XYI;
            Melem(ele,:) = [ne(1), PosN+1, PosN+3]; % sustituyo elemento
            % Este es entre el primer nodo, el punto intermedio entre 1 y 2 y el intermedio entre 1 y3
            
            Melem(PosE+1,:) = [PosN+1, ne(2), PosN+2]; % Este en el punto entre 1 y 2, el nodo 2 y el intermedio entre 2 y 3
            Melem(PosE+2,:) = [PosN+1, PosN+2, PosN+3]; % Este es el central
            Melem(PosE+3,:) = [PosN+3, PosN+2, ne(3)]; %El que queda

            
            % Chequeo de vacios
            PsiABC = N*psie; % Valores de psi para los 3 nodos adicionales
            
            % Elemento sustituido
            if (psie(1) + PsiABC(1) + PsiABC(3))<0
              Vac(ele) = 1;
            end
            % Primer elemento nuevo
            if (PsiABC(1) + psie(2) + PsiABC(2))<0
              Vac(PosE+1)=1;
            end
            % Segundo
            if (PsiABC(1) + PsiABC(2)+PsiABC(3))<0
              Vac(PosE+2)=1;
            end
            % Tercero
            if (PsiABC(3) + PsiABC(2)+psie(3))<0
              Vac(PosE+3)=1;
            end
            
            
            PosE = PosE+3;
            
            % Eta y Xi elemento interfase
            ETI = [etaR1; etaR2; etaR3];
            XII = [xiR1; xiR2; xiR3];
            
            % Valores de las aristas en elemnto interfase
            AR1 = [Ar1R1;Ar1R2;Ar1R3];
            AR2 = [Ar2R1;Ar2R2;Ar2R3];
            
            N = zeros(1,9);
            
            % Valores nodales
            UXX =  [ES.U(2*ne-1); ES.U( reshape( ES.NodAriX(ne,ArX)  , [6,1]) )]; 
            UXY = [ES.U(2*ne); ES.U( reshape( ES.NodAriX(ne,ArX) +1 , [6,1]) )]; 

            for nod = 1:3
                % Coeficientes de la funcion de forma de nodo extendido
                et = ETI(nod); xi = XII(nod);
                N(1) = 1-et-xi; N(2) = et; N(3) = xi;
                
                % Calculo de las g
                g1 = AR1(nod); %g vale simplemente la auxiliar arista en elemento inferfase
                g2 = AR2(nod);
                N(4:6)= N(1:3)*g1;
                N(7:9) = N(1:3)*g2;
                
                UIX = N*UXX;
                UIY = N*UXY;
                
                UiX(PosN+nod) = UIX;
                UiY(PosN+nod) = UIY;
            end
            
            PosN = PosN + 3;
        end
    end
    
    % Melem = Melem(IN>0,:);
    % TR = triangulation(Melem,Mnodo(:,1),Mnodo(:,2),1+0*Mnodo(:,2));
       
   if vacum==0 % En caso de querer borrar los elementos vacios para el plot
     Melem(Vac==1,:)=[]; % Se los elimina de la matriz de elementos para el trisurf
   end
   
    % Escribire en la figura abierta
    
    if bord==0
      Opcion='none';
    else
      Opcion='k';
    end
    
    hold on;
    trisurf(Melem,Mnodo(:,1)+Scale*UiX,Mnodo(:,2)+Scale*UiY,sqrt(UiX.^2+UiY.^2),sqrt(UiX.^2+UiY.^2),'FaceColor','interp','EdgeColor',Opcion);
    
elseif strcmp(tipo,'Psi') % Level Set
  
    Mnodo = [ES.Mnodo(:,2:3); zeros(sum(ES.EI)*3,2)]; % Se inicia matriz de nodos agregando 3 filas por cada elemento extendido
    Melem = [ES.Melem(:,3:5); zeros(sum(ES.EI)*3,3)]; % Idem para la matriz de elementos
    % Como la idea es agregar algunos mas a medida que se encuentren interfases
    % Estos son indices actuales de la posicion ingresada de nodo
    % y de elemento respectivamente
    PosN = ES.Nnodo; 
    PosE = ES.Nelem;
    Ui = [ES.psi; zeros(sum(ES.EI)*3,1)]; % Se inicia vector de la variable a pintar
    
    
    for ele = 1:ES.Nelem
        % Nodos
        ne = ES.Melem(ele,3:5);
        
        % Geometria del elemento
        Xe = ES.Mnodo(ne,2);
        Ye = ES.Mnodo(ne,3);
        XYe=[Xe,Ye];
        
        % Funcion de nivel en el elemento
        psie = ES.psi(ne);
        
        
        
        % ---------------------------------------------------------------------
        % Elemento con interfase ----------------------------------------------
        if ES.EI(ele)
            

            % PUNTOS INTERMEDIOS PARA DIVISION EN SUBELEMENTOS PARA INTEGRACION
            % Y DEFINICION DE LAS "g"

            % Punto entre nodo 1 y 2. En este xi = 0
            if psie(1)*psie(2)<0 
                etaR1 = psie(1)/(psie(1)-psie(2));

            else
                etaR1 = 0.5;
            end
            xiR1=0;

            % Punto entre nodo 2 y 3. En este xi = 1-eta
            if psie(2)*psie(3)<0
                etaR2 = psie(3)/(psie(3)-psie(2));
            else
                etaR2 = 0.5;
            end  
            xiR2=1-etaR2;

            % Punto entre nodo 3 y 1. En este eta=0.        
            if psie(3)*psie(1)<0
                xiR3 = psie(1)/(psie(1)-psie(3));
            else
                xiR3 = 0.5;
            end  
            etaR3=0;
            
            % Nodos adicionales
            N = zeros(3,3);
            N(1,1) = 1-etaR1-xiR1; N(1,2) = etaR1; N(1,3) = xiR1;
            N(2,1) = 1-etaR2-xiR2; N(2,2) = etaR2; N(2,3) = xiR2;
            N(3,1) = 1-etaR3-xiR3; N(3,2) = etaR3; N(3,3) = xiR3;
            XYI = N*XYe; % Coordenadas de los 3 nodos adicionales
            Mnodo(PosN+1:PosN+3,:) = XYI;
            Melem(ele,:) = [ne(1), PosN+1, PosN+3]; % sustituyo elemento
            % Este es entre el primer nodo, el punto intermedio entre 1 y 2 y el intermedio entre 1 y3
            
            Melem(PosE+1,:) = [PosN+1, ne(2), PosN+2]; % Este en el punto entre 1 y 2, el nodo 2 y el intermedio entre 2 y 3
            Melem(PosE+2,:) = [PosN+1, PosN+2, PosN+3]; % Este es el central
            Melem(PosE+3,:) = [PosN+3, PosN+2, ne(3)]; %El que queda

            
            
            PosE = PosE+3;
            
            % Eta y Xi elemento interfase
            ETI = [etaR1; etaR2; etaR3];
            XII = [xiR1; xiR2; xiR3];
            
            N = zeros(1,3);
            
            % Valores nodale
            UX = ES.psi(ne);

            for nod = 1:3
                % Coeficientes de la funcion de forma de nodo extendido
                et = ETI(nod); xi = XII(nod);
                N(1) = 1-et-xi; N(2) = et; N(3) = xi;
              
                
                UI  = N*UX;
                
                Ui(PosN+nod) = UI;
            end
            PosN = PosN + 3;
        end
    end
    

   % Escribire en la figura abierta
    if bord==0
      Opcion='none';
    else
      Opcion='k';
    end
    hold on;
    trisurf(Melem,Mnodo(:,1),Mnodo(:,2),Ui,Ui,'FaceColor','interp','EdgeColor',Opcion);
    
elseif strcmp(tipo,'Mat') % Zonas con material
  
    Mnodo = [ES.Mnodo(:,2:3); zeros(sum(ES.EI)*3,2)]; % Se inicia matriz de nodos agregando 3 filas por cada elemento extendido
    Melem = [ES.Melem(:,3:5); zeros(sum(ES.EI)*3,3)]; % Idem para la matriz de elementos
    % Como la idea es agregar algunos mas a medida que se encuentren interfases
    % Estos son indices actuales de la posicion ingresada de nodo
    % y de elemento respectivamente
    PosN = ES.Nnodo; 
    PosE = ES.Nelem;
     Ui = [ES.psi; zeros(sum(ES.EI)*3,1)]; % Se inicia vector de la variable a pintar
    Vac = zeros(size(Melem,1),1); % Inicializado vector logico para chequeo de vacios
    
    for ele = 1:ES.Nelem
        % Nodos
        ne = ES.Melem(ele,3:5);
        
        % Geometria del elemento
        Xe = ES.Mnodo(ne,2);
        Ye = ES.Mnodo(ne,3);
        XYe=[Xe,Ye];
        
        % Funcion de nivel en el elemento
        psie = ES.psi(ne);
        
		
        if ~ES.EI(ele) && mean(psie)<0 % Si no hay interfaz es facil chequear si es vacio o material
          % Ademas, si era material sin interfaz no hay que hacer nada en Vac,
            Vac(ele) = 1;
          
        end
        
        % ---------------------------------------------------------------------
        % Elemento con interfase ----------------------------------------------
        if ES.EI(ele)
            
            
            % PUNTOS INTERMEDIOS PARA DIVISION EN SUBELEMENTOS PARA INTEGRACION
            % Y DEFINICION DE LAS "g"

            % Punto entre nodo 1 y 2. En este xi = 0
            if psie(1)*psie(2)<0 
                etaR1 = psie(1)/(psie(1)-psie(2));

            else
                etaR1 = 0.5;
            end
            xiR1=0;

            % Punto entre nodo 2 y 3. En este xi = 1-eta
            if psie(2)*psie(3)<0
                etaR2 = psie(3)/(psie(3)-psie(2));
            else
                etaR2 = 0.5;
            end  
            xiR2=1-etaR2;

            % Punto entre nodo 3 y 1. En este eta=0.        
            if psie(3)*psie(1)<0
                xiR3 = psie(1)/(psie(1)-psie(3));
            else
                xiR3 = 0.5;
            end  
            etaR3=0;
            
            
            % Nodos adicionales
            N = zeros(3,3);
            N(1,1) = 1-etaR1-xiR1; N(1,2) = etaR1; N(1,3) = xiR1;
            N(2,1) = 1-etaR2-xiR2; N(2,2) = etaR2; N(2,3) = xiR2;
            N(3,1) = 1-etaR3-xiR3; N(3,2) = etaR3; N(3,3) = xiR3;
            XYI = N*XYe; % Coordenadas de los 3 nodos adicionales
            Mnodo(PosN+1:PosN+3,:) = XYI;
            Melem(ele,:) = [ne(1), PosN+1, PosN+3]; % sustituyo elemento
            % Este es entre el primer nodo, el punto intermedio entre 1 y 2 y el intermedio entre 1 y3
            
            Melem(PosE+1,:) = [PosN+1, ne(2), PosN+2]; % Este en el punto entre 1 y 2, el nodo 2 y el intermedio entre 2 y 3
            Melem(PosE+2,:) = [PosN+1, PosN+2, PosN+3]; % Este es el central
            Melem(PosE+3,:) = [PosN+3, PosN+2, ne(3)]; %El que queda

            
            
            
            % Chequeo de vacios
            PsiABC = N*psie; % Valores de psi para los 3 nodos adicionales
            
            % Elemento sustituido
            if (psie(1) + PsiABC(1) + PsiABC(3))<0
              Vac(ele) = 1;
            end
            % Primer elemento nuevo
            if (PsiABC(1) + psie(2) + PsiABC(2))<0
              Vac(PosE+1)=1;
            end
            % Segundo
            if (PsiABC(1) + PsiABC(2)+PsiABC(3))<0
              Vac(PosE+2)=1;
            end
            % Tercero
            if (PsiABC(3) + PsiABC(2)+psie(3))<0
              Vac(PosE+3)=1;
            end
            
            
            PosE = PosE+3;
             % Eta y Xi elemento interfase
            ETI = [etaR1; etaR2; etaR3];
            XII = [xiR1; xiR2; xiR3];
            
            N = zeros(1,3);
            
            % Valores nodale
            UX = ES.psi(ne);

            for nod = 1:3
                % Coeficientes de la funcion de forma de nodo extendido
                et = ETI(nod); xi = XII(nod);
                N(1) = 1-et-xi; N(2) = et; N(3) = xi;
              
                
                UI  = N*UX;
                
                Ui(PosN+nod) = UI;
            end
            
            PosN = PosN + 3;
        end
    end

   % Escribire en la figura abierta
    if bord==0
      Opcion='none';
    else
      Opcion='k';
    end

    hold on;
    MelemAux=Melem(Vac==0,:);
	
    if size(MelemAux,1)>0
      trisurf(MelemAux,Mnodo(:,1),Mnodo(:,2),Ui,Ui,'FaceColor','k','EdgeColor',Opcion);
    end
    
    MelemAux=Melem(Vac==1,:);
	
    if size(MelemAux,1)>0
       trisurf(MelemAux,Mnodo(:,1),Mnodo(:,2),Ui,Ui,'FaceColor','w','EdgeColor',Opcion);
    end
    
elseif strcmp(tipo,'DTC') % Derivada topologica de la complacencia
  
    Mnodo = ES.Mnodo(:,2:3); % Se inicia matriz de nodos
    Melem = ES.Melem(:,3:5); zeros(sum(ES.EI)*3,3); % Idem para la matriz de elementos
    Ui = ES.DTC; % Se inicia vector de la variable a pintar

   % Escribire en la figura abierta
    if bord==0
      Opcion='none';
    else
      Opcion='k';
    end
    hold on;
    trisurf(Melem,Mnodo(:,1),Mnodo(:,2),Ui,Ui,'FaceColor','interp','EdgeColor',Opcion);
	
elseif strcmp(tipo,'DTT') % Derivada topologica del funcional de tensiones
  
    Mnodo = ES.Mnodo(:,2:3); % Se inicia matriz de nodos
    Melem = ES.Melem(:,3:5); zeros(sum(ES.EI)*3,3); % Idem para la matriz de elementos
    Ui = ES.DTT; % Se inicia vector de la variable a pintar

   % Escribire en la figura abierta
    if bord==0
      Opcion='none';
    else
      Opcion='k';
    end
    hold on;
    trisurf(Melem,Mnodo(:,1),Mnodo(:,2),Ui,Ui,'FaceColor','interp','EdgeColor',Opcion);
    
elseif startsWith(tipo,'Sigma')
      
    Mnodo = [ES.Mnodo(:,2:3); zeros(sum(ES.EI)*3,2)]; % Se inicia matriz de nodos agregando 3 filas por cada elemento extendido
    Melem = [ES.Melem(:,3:5); zeros(sum(ES.EI)*3,3)]; % Idem para la matriz de elementos
    % Como la idea es agregar algunos mas a medida que se encuentren interfases
    % Estos son indices actuales de la posicion ingresada de nodo
    % y de elemento respectivamente
    PosN = ES.Nnodo; 
    PosE = ES.Nelem;
    Vac = zeros(size(Melem,1),1); % Inicializado vector logico para chequeo de vacios

    SumVP = zeros(size(Mnodo,1),1); % Se inicializa un vector de la suma de los VP.
    EleCounter = zeros(size(Mnodo,1),1); % Y un vector que cuenta los elementos por nodo  para calcular la media.
    for ele = 1:ES.Nelem
        % Nodos
        ne = ES.Melem(ele,3:5);
        
        % Geometria del elemento
        Xe = ES.Mnodo(ne,2);
        Ye = ES.Mnodo(ne,3);
        XYe=[Xe,Ye];
        
        % Funcion de nivel en el elemento
        psie = ES.psi(ne);
        
		  
        % Derivadas de las funciones de forma segun coord. intrinsecas
        dN_detachi=[-1 1 0; -1 0 1];

        % Matriz Jacobiana
        J = dN_detachi*XYe;

        % Derivadas de las funciones de forma segun x e y.
        dN_dxy = J\dN_detachi; % Se hace con J\ en vez de lo usual de inv(J)* pues 
        % Matlab dice que asi es mas rapido.
    
    
        if ~ES.EI(ele) && mean(psie)<0 % Si no hay interfaz es facil chequear si es vacio o material
          % Ademas, si era material sin interfaz no hay que hacer nada en Vac,
            Vac(ele) = 1;
          
        end
        
        % ---------------------------------------------------------------------
        % Elemento con interfase ----------------------------------------------
        if ES.EI(ele)
            
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
            xiR1=0;

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
            xiR2=1-etaR2;

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
            etaR3=0;

            % Nodos adicionales
            N = zeros(3,3);
            N(1,1) = 1-etaR1-xiR1; N(1,2) = etaR1; N(1,3) = xiR1;
            N(2,1) = 1-etaR2-xiR2; N(2,2) = etaR2; N(2,3) = xiR2;
            N(3,1) = 1-etaR3-xiR3; N(3,2) = etaR3; N(3,3) = xiR3;
            XYI = N*XYe; % Coordenadas de los 3 nodos adicionales
            Mnodo(PosN+1:PosN+3,:) = XYI;
            Melem(ele,:) = [ne(1), PosN+1, PosN+3]; % sustituyo elemento
            % Este es entre el primer nodo, el punto intermedio entre 1 y 2 y el intermedio entre 1 y3
            
            Melem(PosE+1,:) = [PosN+1, ne(2), PosN+2]; % Este en el punto entre 1 y 2, el nodo 2 y el intermedio entre 2 y 3
            Melem(PosE+2,:) = [PosN+1, PosN+2, PosN+3]; % Este es el central
            Melem(PosE+3,:) = [PosN+3, PosN+2, ne(3)]; %El que queda
            
            % Chequeo de vacios
            PsiABC = N*psie; % Valores de psi para los 3 nodos adicionales
            
             
            % Elemento sustituido
            if (psie(1) + PsiABC(1) + PsiABC(3))<0
              Vac(ele) = 1;
            end
            % Primer elemento nuevo
            if (PsiABC(1) + psie(2) + PsiABC(2))<0
              Vac(PosE+1)=1;
            end
            % Segundo
            if (PsiABC(1) + PsiABC(2)+PsiABC(3))<0
              Vac(PosE+2)=1;
            end
            % Tercero
            if (PsiABC(3) + PsiABC(2)+psie(3))<0
              Vac(PosE+3)=1;
            end
            
            % Las matrices Ar1SE y Ar2SE son distintas al resto de los
            % codigos porque son acordes a los subelementos descritos
            % arriba.
            
             Ar1SE = [0      Ar1R1      Ar1R1     Ar1R2;
                     Ar1R1  0   Ar1R2 0;
                     Ar1R3  Ar1R2  Ar1R3 Ar1R3 ]; % Formato similar a los anteriores para los valores de g1 en los vertices de los subtriangulos
    
            Ar2SE = [0      Ar2R1      Ar2R1     Ar2R2;
                     Ar2R1  0      Ar2R2 0;
                     Ar2R3  Ar2R2  Ar2R3 Ar2R3 ]; % Formato similar a los anteriores para los valores de g2 en los vertices de los subtriangulos
             
            
            PosE = PosE+3;
    
            etaG = [1/3 0.6 0.2 0.2]; % Coordenadas eta de los puntos.
            xiG = [1/3 0.2 0.6 0.2]; % Coordenadas xi de los puntos
            PesoG = [-0.28125 0.260416666666666667 0.260416666666666667  0.260416666666666667]; % Peso de los puntos.
        
            % Primer nuevo elemento
            etase = [0; etaR1; 0]; %etas de este elemento
            xise = [0; 0; xiR3]; % xis.
            
            Ar1 = Ar1SE(:,1); % Valores de Ar para g1
            Ar2 = Ar2SE(:,1); %Valores de Ar para g2
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
            else
                E=ES.PropMat(2,2); %Young
                nu=ES.PropMat(2,3); %Poisson
            end
            

            if ES.PEL==2     
               E = E/(1-nu^2);
               nu = nu / (1-nu);
            end        
        
            C = (E/(1-nu^2))* [ 1 nu 0;
                            nu 1 0;
                            0 0 0.5*(1-nu)]; % Matriz tensor elastico
              
            % Integracion en puntos de Gauss
            
            for puntoG=1:4 % Para los 3 puntos de Gauss anteriormente definidos
                
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
                 %En realidad la componente del tensor fuera de la diagonal es la
                 %mitad de la asignada en ese vector a la tercer componente
                % pero como luego se va a hacer un producto escalar con un tensor simetrico
                % podemos hacer este producto una vez usando este valor.      
                % Ademas C se arma considerando que se esta considerando la
                % distorsion angular y no la componente del tensor. 
                
                 % Vector de tensiones elementales            
                Sigma=C*Eps;

                if strcmp(tipo, 'SigmaX')
                    Sigma=Sigma(1);
                elseif strcmp(tipo, 'SigmaY')
                    Sigma=Sigma(2);
                elseif strcmp(tipo,'SigmaXY')
                    Sigma=Sigma(3);
                elseif strcmp(tipo,'SigmaMises')
                    Sigma=sqrt( (Sigma(1) + Sigma(2))^2 -3*(Sigma(1)*Sigma(2) -Sigma(3)^2));
                end

                ValPerNodo = ( Nse * Sigma )'; % Valores por nodo en este punto

                SumVP(ne(1)) = SumVP(ne(1)) + 6 * Jse * PesoG(puntoG) * ValPerNodo(1); % Se suma este punto de Gauss
                SumVP(PosN+1) = SumVP(PosN+1) + 6 * Jse * PesoG(puntoG) * ValPerNodo(2);
                SumVP(PosN+3) =  SumVP(PosN+3) + 6 * Jse * PesoG(puntoG) * ValPerNodo(3);
             
                
            end %En for en los puntos de Gauss
            EleCounter(ne(1))=EleCounter(ne(1))+Jse;
            EleCounter(PosN+1)=EleCounter(PosN+1)+Jse;
            EleCounter(PosN+3)=EleCounter(PosN+3)+Jse;
            
            
            
            % Segundo nuevo elemento
            etase = [etaR1; 1; etaR2]; %etas de este elemento
            xise = [0; 0; 1-etaR2]; % xis.
              
            Ar1 = Ar1SE(:,2); % Valores de Ar para g1
            Ar2 = Ar2SE(:,2); %Valores de Ar para g2
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
            else
                E=ES.PropMat(2,2); %Young
                nu=ES.PropMat(2,3); %Poisson
            end
            

            if ES.PEL==2     
               E = E/(1-nu^2);
               nu = nu / (1-nu);
            end        
        
        
            C = (E/(1-nu^2))* [ 1 nu 0;
                            nu 1 0;
                            0 0 0.5*(1-nu)]; % Matriz tensor elastico
              
            % Integracion en puntos de Gauss
            
            for puntoG=1:4 % Para los 3 puntos de Gauss anteriormente definidos
                
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
                 %En realidad la componente del tensor fuera de la diagonal es la
                 %mitad de la asignada en ese vector a la tercer componente
                % pero como luego se va a hacer un producto escalar con un tensor simetrico
                % podemos hacer este producto una vez usando este valor.      
                % Ademas C se arma considerando que se esta considerando la
                % distorsion angular y no la componente del tensor. 
                
                 % Vector de tensiones elementales            
                Sigma=C*Eps;

                if strcmp(tipo, 'SigmaX')
                    Sigma=Sigma(1);
                elseif strcmp(tipo, 'SigmaY')
                    Sigma=Sigma(2);
                elseif strcmp(tipo,'SigmaXY')
                    Sigma=Sigma(3);
                elseif strcmp(tipo,'SigmaMises')
                    Sigma=sqrt( (Sigma(1) + Sigma(2))^2 -3*(Sigma(1)*Sigma(2) -Sigma(3)^2));
                end

                ValPerNodo = ( Nse * Sigma )'; % Valores por nodo en este punto

                SumVP(PosN+1) = SumVP(PosN+1)+ 6 * Jse * PesoG(puntoG) * ValPerNodo(1); % Se suma este punto de Gauss
                SumVP(ne(2)) =  SumVP(ne(2))+ 6 * Jse * PesoG(puntoG) * ValPerNodo(2);
                SumVP(PosN+2) = SumVP(PosN+2)+  6 * Jse * PesoG(puntoG) * ValPerNodo(3);
             
                
            end %En for en los puntos de Gauss
            EleCounter(PosN+1)=EleCounter(PosN+1)+Jse;
            EleCounter(ne(2))=EleCounter(ne(2))+Jse;
            EleCounter(PosN+2)=EleCounter(PosN+2)+Jse;
            
            
            
            % Tercer nuevo elemento
            etase = [etaR1; etaR2; 0]; %etas de este elemento
            xise = [0; 1-etaR2; xiR3]; % xis.
           
            Ar1 = Ar1SE(:,3); % Valores de Ar para g1
            Ar2 = Ar2SE(:,3); %Valores de Ar para g2
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
            else
                E=ES.PropMat(2,2); %Young
                nu=ES.PropMat(2,3); %Poisson
            end
            

            if ES.PEL==2     
               E = E/(1-nu^2);
               nu = nu / (1-nu);
            end        
        
        
            C = (E/(1-nu^2))* [ 1 nu 0;
                            nu 1 0;
                            0 0 0.5*(1-nu)]; % Matriz tensor elastico
              
            % Integracion en puntos de Gauss
            
            for puntoG=1:4 % Para los 3 puntos de Gauss anteriormente definidos
                
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
                 %En realidad la componente del tensor fuera de la diagonal es la
                 %mitad de la asignada en ese vector a la tercer componente
                % pero como luego se va a hacer un producto escalar con un tensor simetrico
                % podemos hacer este producto una vez usando este valor.      
                % Ademas C se arma considerando que se esta considerando la
                % distorsion angular y no la componente del tensor. 
                
                 % Vector de tensiones elementales            
                Sigma=C*Eps;

                if strcmp(tipo, 'SigmaX')
                    Sigma=Sigma(1);
                elseif strcmp(tipo, 'SigmaY')
                    Sigma=Sigma(2);
                elseif strcmp(tipo,'SigmaXY')
                    Sigma=Sigma(3);
                elseif strcmp(tipo,'SigmaMises')
                    Sigma=sqrt( (Sigma(1) + Sigma(2))^2 -3*(Sigma(1)*Sigma(2) -Sigma(3)^2));
                end

                ValPerNodo = ( Nse * Sigma )'; % Valores por nodo en este punto

                SumVP(PosN+1) = SumVP(PosN+1)+ 6 * Jse * PesoG(puntoG) * ValPerNodo(1); % Se suma este punto de Gauss
                SumVP(PosN+2) =  SumVP(PosN+2)+ 6 * Jse * PesoG(puntoG) * ValPerNodo(2);
                SumVP(PosN+3) = SumVP(PosN+3)+  6 * Jse * PesoG(puntoG) * ValPerNodo(3);
             
                
            end %En for en los puntos de Gauss
            EleCounter(PosN+1)=EleCounter(PosN+1)+Jse;
            EleCounter(PosN+2)=EleCounter(PosN+2)+Jse;
            EleCounter(PosN+3)=EleCounter(PosN+3)+Jse;
            
            
            
            % Cuarto (y ultimo) nuevo elemento
            etase = [etaR2; 0; 0]; %etas de este elemento
            xise = [1-etaR2; 1; xiR3]; % xis.
            
            Ar1 = Ar1SE(:,4); % Valores de Ar para g1
            Ar2 = Ar2SE(:,4); %Valores de Ar para g2
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
            else
                E=ES.PropMat(2,2); %Young
                nu=ES.PropMat(2,3); %Poisson
            end
            

            if ES.PEL==2     
               E = E/(1-nu^2);
               nu = nu / (1-nu);
            end        
        
            C = (E/(1-nu^2))* [ 1 nu 0;
                            nu 1 0;
                            0 0 0.5*(1-nu)]; % Matriz tensor elastico
              
            % Integracion en puntos de Gauss
            
            for puntoG=1:4 % Para los 3 puntos de Gauss anteriormente definidos
                
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
                 %En realidad la componente del tensor fuera de la diagonal es la
                 %mitad de la asignada en ese vector a la tercer componente
                % pero como luego se va a hacer un producto escalar con un tensor simetrico
                % podemos hacer este producto una vez usando este valor.      
                % Ademas C se arma considerando que se esta considerando la
                % distorsion angular y no la componente del tensor. 
                
                 % Vector de tensiones elementales            
                Sigma=C*Eps;

                if strcmp(tipo, 'SigmaX')
                    Sigma=Sigma(1);
                elseif strcmp(tipo, 'SigmaY')
                    Sigma=Sigma(2);
                elseif strcmp(tipo,'SigmaXY')
                    Sigma=Sigma(3);
                elseif strcmp(tipo,'SigmaMises')
                    Sigma=sqrt( (Sigma(1) + Sigma(2))^2 -3*(Sigma(1)*Sigma(2) -Sigma(3)^2));
                end
                

                ValPerNodo = ( Nse * Sigma )'; % Valores por nodo en este punto

                SumVP(PosN+2) = SumVP(PosN+2)+ 6 * Jse * PesoG(puntoG) * ValPerNodo(1); % Se suma este punto de Gauss
                SumVP(ne(3)) =  SumVP(ne(3))+ 6 * Jse * PesoG(puntoG) * ValPerNodo(2);
                SumVP(PosN+3) = SumVP(PosN+3)+  6 * Jse * PesoG(puntoG) * ValPerNodo(3);
             
                
            end %En for en los puntos de Gauss
            EleCounter(PosN+2)=EleCounter(PosN+2)+Jse;
            EleCounter(ne(3))=EleCounter(ne(3))+Jse;
            EleCounter(PosN+3)=EleCounter(PosN+3)+Jse;
            
            
            
            
            
            PosN = PosN + 3;
        
     % --------------------------------------------
    % Elemento vacio -----------------------------
    elseif ES.VA(ele)
        % Asignacion de propiedades del material estructural
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
        
        
                if strcmp(tipo, 'SigmaX')
                    Sigma=Sigma(1);
                elseif strcmp(tipo, 'SigmaY')
                    Sigma=Sigma(2);
                elseif strcmp(tipo,'SigmaXY')
                    Sigma=Sigma(3);
                elseif strcmp(tipo,'SigmaMises')
                    Sigma=sqrt( (Sigma(1) + Sigma(2))^2 -3*(Sigma(1)*Sigma(2) -Sigma(3)^2));
                end
        
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
            

                ValPerNodo = ( N * Sigma )'; % Valores por nodo en este punto
            
            SumVP(ne) =SumVP(ne)+ 6 * PesoG(puntoG) * ValPerNodo; % Se suma este punto de Gauss
            
        end % Fin de integrar en puntos de gauss
        EleCounter(ne)=EleCounter(ne)+ones(3,1);
        
        
        % Material Estructural
        % --------------------
        else
            
            
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
        
        
                if strcmp(tipo, 'SigmaX')
                    Sigma=Sigma(1);
                elseif strcmp(tipo, 'SigmaY')
                    Sigma=Sigma(2);
                elseif strcmp(tipo,'SigmaXY')
                    Sigma=Sigma(3);
                elseif strcmp(tipo,'SigmaMises')
                    Sigma=sqrt( (Sigma(1) + Sigma(2))^2 -3*(Sigma(1)*Sigma(2) -Sigma(3)^2));
                end
        
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
            

                ValPerNodo = ( N * Sigma )'; % Valores por nodo en este punto
            
            SumVP(ne) =SumVP(ne)+ 6 * PesoG(puntoG) * ValPerNodo; % Se suma este punto de Gauss
            
        end % Fin de integrar en puntos de gauss
        EleCounter(ne)=EleCounter(ne)+ones(3,1);
            
            
        end %Fin de tipo de elemento
        
    end
    
    % Resta uniformizar porque hay algunos nodos extendidos que aparecen
    % repetido
    % Eso no es problema en el resto de los tipos pero en tensiones SI LO
    % ES.
    
    PosN = ES.Nnodo+1;
    while PosN<length(SumVP)
        d=zeros(length(SumVP)-PosN,1);
        for i=1:length(d)
            d(i) = (Mnodo(i+PosN,1)-Mnodo(PosN,1))^2 + (Mnodo(i+PosN,2)-Mnodo(PosN,2))^2;
        end
        DeltaDist=0.5*ES.Lx/ES.Nx + 0.5*ES.Ly/ES.Ny;
        indAux=d<DeltaDist*1e-12;
        ind=1:length(d); ind=ind+PosN;
        ind=ind(indAux);
        
        for i=1:length(ind)
            indice=ind(i);
            SumVP(PosN)=SumVP(PosN)+SumVP(indice);
            EleCounter(PosN)=EleCounter(PosN)+EleCounter(indice);
            Melem(Melem==indice)=PosN;
            Melem(Melem>indice)=Melem(Melem>indice)-1;
        end
        
        
        SumVP(ind)=[];
        EleCounter(ind)=[];
        Mnodo(ind,:)=[];
        PosN=PosN+1;
    end
    
    
    Ui = SumVP./EleCounter; % Se calcula la media 
    
    
   if vacum==0 % En caso de querer borrar los elementos vacios para el plot
     Melem(Vac==1,:)=[]; % Se los elimina de la matriz de elementos para el trisurf
   end
   

   % Escribire en la figura abierta
    if bord==0
      Opcion='none';
    else
      Opcion='k';
    end
    hold on;
    trisurf(Melem,Mnodo(:,1),Mnodo(:,2),Ui,Ui,'FaceColor','interp','EdgeColor',Opcion);
elseif strcmp(tipo,'B')
    
% Eta y XI en los nodos del triangulo
ETN = [0 1 0];
XIN = [0 0 1];
    hold on;
    
    for ele = 1:ES.Nelem
        
        % ---------------------------------------------------------------------
        % Elemento con interfase ----------------------------------------------
        if ES.EI(ele)
            
            % Nodos
            ne = ES.Melem(ele,3:5);
            
            % Funcion de nivel
            phe = ES.psi(ne);
            
            % Geometria elemento
            XYe = ES.Mnodo(ne,2:3);
            
            pc = 0;
            etc = zeros(1,2);
            xic = zeros(1,2);
            
            % Posicion de interfase
            if (phe(1)*phe(2)<0)
                % corta interfase
                pc = pc + 1;
                dp = phe(2)-phe(1);
                etc(pc) = phe(2)/dp*ETN(1)-phe(1)/dp*ETN(2);
                xic(pc) = phe(2)/dp*XIN(1)-phe(1)/dp*XIN(2);
            end
            if (phe(2)*phe(3)<0)
                % corta interfase
                pc = pc + 1;
                dp = phe(3)-phe(2);
                etc(pc) = phe(3)/dp*ETN(2)-phe(2)/dp*ETN(3);
                xic(pc) = phe(3)/dp*XIN(2)-phe(2)/dp*XIN(3);
            end
            if (phe(3)*phe(1)<0)
                % corta interfase
                pc = pc + 1;
                dp = phe(1)-phe(3);
                etc(pc) = phe(1)/dp*ETN(3)-phe(3)/dp*ETN(1);
                xic(pc) = phe(1)/dp*XIN(3)-phe(3)/dp*XIN(1);
            end
            
            
            % Nodos adicionales
            N = zeros(2,3);
            N(:,1) = 1-etc-xic; N(:,2) = etc; N(:,3) = xic;
            XYI = N*XYe;
            
            % plot3(XYI(:,1), XYI(:,2),  15+0*XYI(:,2), 'r', 'linewidth',2)
            plot(XYI(:,1), XYI(:,2), 'r', 'linewidth',2)
        end
    end    
    
end


end