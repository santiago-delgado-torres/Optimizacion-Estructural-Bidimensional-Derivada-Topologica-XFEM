function ES = DirichletXFEM(ES)

% Rutina para incorporar los efectos de los desplazamientos prescritos
%
% Devuelve o actualiza:
%
% Ktotal: Matriz del sistema (Puede ser rotada por los inclinados
% Fext: Vector de fuerzas externas (Puede ser rotado por los inclinados y modificado por
%       prescripciones de desplazamientos no nulas)
% gdlfij: Grados de libertad prescritos. Vector logico de igual cantidad de entradas que GL
% Upres: Valores para los grados de libertad prescritos. Es un vector de entradas nulas
%        de igual cantidad de entradas que GL la estructura pero en los grados de libertad prescrito
%        posee estos mismos.
% GdLiRot: Vector de la mitad de entradas que los grados de libertad. Muestra el angulo
%          de rotacion actual de las variables del nodo (tradicional o no).
%          El angulo es positivo antihorario desde el eje x.

% =========================================================================
% === Inicializacion de Variables =========================================
% =========================================================================

ES.gdlfij = zeros(2*ES.Nnodo + ES.NGLX , 1 ); 
ES.Upres = zeros(2*ES.Nnodo + ES.NGLX , 1);
ES.GdLiRot = zeros(ES.Nnodo + ES.NGLX/2, 1 ); % Afortunadamente ES.NGLX siempre es par 

% =========================================================================
% === Apoyos Puntuales Fijos ==============================================
% =========================================================================


for i = 1:size(ES.CB.Dir.Nod.Fijo,1) % Recorremos todas las entradas de la matriz asociada
   
   % Estos no requieren ningun tratamiento especial mas que registrar los valores
   % Se recuerda que estos NO afectana los grados de libertad extendidos
   
   ES.gdlfij( (2*ES.CB.Dir.Nod.Fijo(i,1)-1):(2*ES.CB.Dir.Nod.Fijo(i,1)) ) = ones(2,1);  %Se fijan estos 2
   
   ES.Upres(2*ES.CB.Dir.Nod.Fijo(i,1)-1) = ES.CB.Dir.Nod.Fijo(i,2); % Se guarda el de X
   ES.Upres(2*ES.CB.Dir.Nod.Fijo(i,1)  ) = ES.CB.Dir.Nod.Fijo(i,3); % Se guarda el Y   
   
end

% ==========================================================================
% === Apoyos Puntuales Deslizantes =========================================
% ==========================================================================

for i = 1:size(ES.CB.Dir.Nod.Desl,1) % Se recorren todas las entradas de la matriz asociada
      
   % Notese que al ser puntuales no afectan a los Grados de libertad extendidos
   
   nod = ES.CB.Dir.Nod.Desl(i,1) ;  % Numero de nodo afectado
   theta = ES.CB.Dir.Nod.Desl(i,2) - ES.GdLiRot( nod ); % Angulo de los nuevos ejes respecto al anterior 
                                                        % El e_1 de la nueva base es en el cual desliza.
                                                        
   if theta==0 || abs(theta)==180 %Caso en el cual desliza segun el eje e_1
       % Modificacion de Upres
       % ---------------------
       if theta==0
           ES.Upres( 2*nod ) = ES.CB.Dir.Nod.Desl(i,3);
       else
           ES.Upres( 2*nod ) = - ES.CB.Dir.Nod.Desl(i,3); % Si fuera ± 180 el e2 prescrito es -e2, entonces el signo
       end
       % Modificacion de gdlfig
       % ----------------------
       ES.gdlfij(2*nod) = 1;
   elseif abs(theta)==90 || abs(theta)==270
       % Modificacion de Upres
       % ---------------------
       if theta==90 || theta==-270
           ES.Upres( 2*nod - 1 ) = - ES.CB.Dir.Nod.Desl(i,3); % En estos angulos el prescrito es -e1 actual
       else
           ES.Upres( 2*nod - 1) =  ES.CB.Dir.Nod.Desl(i,3); % En estos angulos el prescrito es e1 actual
       end       
       % Modificacion de gdlfig
       % ----------------------
       ES.gdlfij( 2*nod -1 ) = 1;
   else
       % Este lleva un trabajo levemente mayor.
       % Al tener apoyos inclinados (al menos el caso general). Lo que suele hacerse
       % es rotar las 2 filas y las 2 columnas de la matriz asociadas al desplazamiento
       % del nodo con apoyo inclinado y luego que el vector de fuerzas externas y el 
       % vector de desplazamientos esten expresados en una base solidaria a la inclinacion del apoyo


       Tr = [ cosd(theta) , -sind(theta) ;
              sind(theta) , cosd(theta) ]; %Matriz de rotacion

       Tr = sparse(Tr) ; %Pasada a dispersa para que Matlab/Octave pueda hacer la cuenta de rotacion de la matriz

       % Modificacion de Upres
       %-----------------------

       % Este paso importa porque si por algun otro contexto se prescribio el desplazamiento que con este queda deslizante
       % no se quiere perder esta prescripcion original

       ES.Upres( (2*nod) + (-1:0) ) = Tr' * ES.Upres( (2*nod) + (-1:0) ); %Tr pasa de nuevos a anteriores
                                                                          % Tr' pasa de anteriores a nuevos

       % Se presciribe el que esta prescripto:
       ES.Upres( 2*nod ) = ES.CB.Dir.Nod.Desl(i,3);

       % Modificacion de gdlfig
       % ----------------------

       % Se distinguen dos casos.
       % 
       % 1) Que este nodo ya tuviera alguna parte prescrita
       %    
       %    En este caso, a menos que el prescrito coincida con el que prescribe este apoyo
       %    Que seria raro en esta seccion pero ya sirve para los de borde 
       %    Queda prescrito TODO el nodo
       %    Si juuuusto el prescrito coincide con el que es prescrito
       %    Entonces aplica lo mismo que el segundo caso.
       % 
       % 2) Que el nodo no tuviera nada prescrito
       % 
       %    En este caso solamente se debe prescribir el 2*nod

       if ES.gdlfij( 2*nod - 1).^2 + ES.gdlfij( 2*nod ).^2 ~= 0 % Deben ser los dos nulos (no prescritos) para no entrar en este if
          if ES.gdlfij( 2*nod - 1 )*ES.gdlfij( 2*nod) ==0 % Si alguno de los dos no estaba prescrito
            if ES.gdlfij(2*nod - 1)==1 % El primero estaba prescrito
              if theta==90 || theta==270 % Si el primero ahora es el segundo (y este es el que queremos prescribir
                ES.gdlfij(2*nod) = 1;
                ES.gdlfij(2*nod - 1 )=0;
                % Como el 1 original pasa a ser el 2
                % Lo dejamos prescrito prescribiendo el 2
                %
                % Como el 2 original (libre) pasa a ser el 1, y este apoyo no prescribe el 1
                % Desprescribimos el 1.
              else
                % Queda prescrito todo.
                % El 1 ya esta
                % El 2:
                ES.gdlfij(2*nod) = 1;
              end
            else % El que estaba prescrito era el 2
              if theta~=0 && theta~=180 % A menos que los ejes no cambien, solo evenutlamente de signo
                % Hay que prescribir todo
                % el 2 ya esta
                ES.gdlfij(2*nod-1) = 1;
              end %No aplica else. Si el 2 estaba restringido, va a seguir y el 1 no lo estaba y no debe seguir
            end
          end % No aplica un else pues si ya estaban prescritos... que sigan
       else 
          ES.gdlfij( 2*nod ) = 1; % Se prescribe el segundo grado de libertad solamente 
       end  


       % Modificacion de Fext
       %---------------------

       ES.Fext( (2*nod) + (-1:0) ) = Tr' * ES.Fext( (2*nod) + (-1:0) ); % Idem explicacion

       % Modificacion de Ktotal
       %------------------------

       Kaux = ES.Ktotal ; % Se la define como una nueva variable.
       % Esto es para no andar modificando mucho en la variable de la struct que suele ser mas lento

       % Rotacion de la columna
       for j = 1 : (ES.Nnodo +  ES.NGLX/2) % Como se rotan de a matrices 2x2, el numero de matrices a rotar es la mitad de las filas

         Kaux( (j*2-1):(j*2) , 2*nod + (-1:0)  ) = Kaux( (j*2-1):(j*2) , 2*nod + (-1:0) ) * Tr ;

       end

       % Rotacion de la fila
       for j = 1 : (ES.Nnodo +  ES.NGLX/2) % Como se rotan de a matrices 2x2, el numero de matrices a rotar es la mitad de las filas

         Kaux( 2*nod + (-1:0) , (j*2-1):(j*2) ) = Tr' * Kaux( 2*nod + (-1:0) , (j*2-1):(j*2) );

       end   

       % Actualizacion de K

       ES.Ktotal = Kaux;  

       % Actualizacion de GdLiRot 
       % ------------------------

       ES.GdLiRot( nod ) = ES.CB.Dir.Nod.Desl(i,2);
   end
  
end % Se termina de recorrer la matriz


% =========================================================================
% === Apoyos Lineales Fijos ==============================================
% =========================================================================

for i = 1:size(ES.CB.Dir.Lin.Fijo,1) % Recorremos todas las entradas de la matriz asociada
   
   % Estos no requieren ningun tratamiento especial mas que registrar los valores
   % Estos si pueden afectar a los grados de libertd extendidos (Los dejan nulos si y solo si hay interfaz en el borde)
   
    ele = ES.CB.Dir.Lin.Fijo(i,1) ; % Variable auxiliar con el numero de elemento para no tener que volver a escribir todo
    
    ne  = ES.Melem(ele,[3:5,3]); % Los nodos del elemento. Lo de repetir el 3 es para un truquito con el siguiente paso
    
    nb = ne( ES.CB.Dir.Lin.Fijo(i,2):(ES.CB.Dir.Lin.Fijo(i,2)+1) ); %Nodos del borde apoyado
    
    SeRestringeGX = ( (ES.psi( nb(1) ) * ES.psi( nb(2) )) < 0 ) ; 
    % Variable logica, verifica que haya intefase en el borde.
    % Si no la hay no se van a tocar los extendidos pero si la hay entonces si
           
    % GL tradicionales:
    % -----------------
    
    ES.gdlfij( 2*nb - 1 ) = ones(2,1); %Se fijan los tradicionales segun x
    ES.gdlfij( 2*nb     ) = ones(2,1); % Idem para los segun y
    
    % Como es posible que algun GL tradicional este rotado (aunque seria raro) vamos a rotar las prescripciones
    
    for nod = nb % se recorren los 2 nodos
      theta = ES.GdLiRot( nod ) ; %Angulo del nodo
      
      Tr = [ cosd(theta) , -sind(theta) ;
             sind(theta) , cosd(theta) ]; %Matriz de rotacion

      Prex =  ES.CB.Dir.Lin.Fijo(i,3:4)' ; % VAlores XY a prescribir

      Prex = Tr'*Prex; % Se rotan a la base del nodo

      ES.Upres( 2*nod+(-1:0) ) = Prex; % Se asignan esos valores al nodo 
  
    end     
    
    
    % GL extendidos:
    %----------------
    if SeRestringeGX
      
      % De manera análoga a las fuerzas, solamente importa los grados
      % extendidos del que tiene la arista en el borde
      
      ArX = ES.EGLX(ele,:); %Numero de aristas extendidas de este elemento
      NArX = ES.AriX(ArX,:); %Numero de nodo de esas aristas
      if ( nb(1)==NArX(1,1) || nb(1)==NArX(1,2) ) && ( nb(2)==NArX(1,1) || nb(2) == NArX(1,2) ) % La fuerza es en la arista extendida 1 
          AriX = 1; 
      else % Es la arista 2
          AriX=2;
      end
        
      nbx = zeros(4,1); %Se inicializa lista con lo GL extendidos a restringir
      
      nbx(1:2) = ES.NodAriX(nb(1),ArX(AriX))' + (0:1); % Los del primer nodo
      
      nbx(3:4) = ES.NodAriX(nb(2),ArX(AriX))'+ (0:1); % Los del segundo nodo
      
      ES.gdlfij( nbx ) = ones(4,1); % Se fijan estos 4
      
      % Aunque estuvieran rotados se los prescribe como nulos
      % Asi que no hay que hacer mucho mas que:
      
      ES.Upres( nbx ) = zeros(4,1); 
      
    end
   
end

% ==========================================================================
% === Apoyos Lineales Deslizantes ==========================================
% ==========================================================================

for i = 1:size(ES.CB.Dir.Lin.Desl,1) % Se recorren todas las entradas de la matriz asociada
  
   
   % Estos no requieren ningun tratamiento especial mas que registrar los valores
   % Estos si pueden afectar a los grados de libertd extendidos (Los dejan nulos si y solo si hay interfaz en el borde)
   
   ele = ES.CB.Dir.Lin.Desl(i,1) ; % Variable auxiliar con el numero de elemento para no tener que volver a escribir todo
    
   ne  = ES.Melem(ele,[3:5,3]); % Los nodos del elemento. Lo de repetir el 3 es para un truquito con el siguiente paso
    
   nb = ne( ES.CB.Dir.Lin.Desl(i,2):(ES.CB.Dir.Lin.Desl(i,2)+1) ); %Nodos del borde apoyado
    
   SeRestringeGX = ( (ES.psi( nb(1) ) * ES.psi( nb(2) )) < 0 ) ; 
   % Variable logica, verifica que haya intefase en el borde.
   % Si no la hay no se van a tocar los extendidos pero si la hay entonces si
   
   % Solamente interesa ademas los extendidos de la arista que coincida con
   % el borde:
   if SeRestringeGX
      ArX = ES.EGLX(ele,:); %Numero de aristas extendidas de este elemento
      NArX = ES.AriX(ArX,:); %Numero de nodo de esas aristas
      if ( nb(1)==NArX(1,1) || nb(1)==NArX(1,2) ) && ( nb(2)==NArX(1,1) || nb(2) == NArX(1,2) ) % La fuerza es en la arista extendida 1 
          AriX = 1; 
      else % Es la arista 2
          AriX=2;
      end
   end
    
   % En este caso va a ser mas facil hacer un for en los dos nodos
   
   for nod = nb
      
      % Primero GL tradicional
      theta = ES.CB.Dir.Lin.Desl(i,3) - ES.GdLiRot( nod ); % Angulo de los nuevos ejes respecto al anterior 
                                                        % El e_1 de la nueva base es en el cual desliza.
                                                        
       if theta==0 || abs(theta)==180 %Caso en el cual desliza segun el eje e_1
           % Modificacion de Upres
           % ---------------------
           if theta==0
               ES.Upres( 2*nod ) = ES.CB.Dir.Lin.Desl(i,4);
           else
               ES.Upres( 2*nod ) = - ES.CB.Dir.Lin.Desl(i,4); % Si fuera ± 180 el e2 prescrito es -e2, entonces el signo
           end
           % Modificacion de gdlfig
           % ----------------------
           ES.gdlfij(2*nod) = 1;
       elseif abs(theta)==90 || abs(theta)==270
           % Modificacion de Upres
           % ---------------------
           if theta==90 || theta==-270
               ES.Upres( 2*nod - 1 ) = - ES.CB.Dir.Lin.Desl(i,4); % En estos angulos el prescrito es -e1 actual
           else
               ES.Upres( 2*nod - 1) =  ES.CB.Dir.Lin.Desl(i,4); % En estos angulos el prescrito es e1 actual
           end       
           % Modificacion de gdlfig
           % ----------------------
           ES.gdlfij( 2*nod -1 ) = 1;
       else


          Tr = [ cosd(theta) , -sind(theta) ;
              sind(theta) , cosd(theta) ]; %Matriz de rotacion

          Tr = sparse(Tr) ; %Pasada a dispersa para que Matlab/Octave pueda hacer la cuenta de rotacion de la matriz     


          %Modificacion de Upres tradicional

          % Este paso importa pues si por algun otro contexto se prescribio el
          % desplazamiento que con este queda deslizante; no se quiere perder esta
          % prescripcion original

          ES.Upres( (2*nod) + (-1:0) ) = Tr' * ES.Upres( (2*nod) + (-1:0) );

          % Se prescribe el que esta prescripto
          ES.Upres( 2*nod ) = ES.CB.Dir.Lin.Desl(i,4);

          % Modificacion de gdlfij en tradicional
          % Ver justificacion del if en deslizantes nodales

          if ES.gdlfij( 2*nod - 1).^2 + ES.gdlfij( 2*nod ).^2 ~= 0 % Deben ser los dos nulos (no prescritos) para no entrar en este if
              if ES.gdlfij( 2*nod - 1 )*ES.gdlfij( 2*nod) ==0 % Si alguno de los dos no estaba prescrito
                if ES.gdlfij(2*nod - 1)==1 % El primero estaba prescrito
                  if theta==90 || theta==270 % Si el primero ahora es el segundo (y este es el que queremos prescribir
                   ES.gdlfij(2*nod) = 1;
                   ES.gdlfij(2*nod - 1 )=0;
                   % Como el 1 original pasa a ser el 2
                   % Lo dejamos prescrito prescribiendo el 2
                    %
                   % Como el 2 original (libre) pasa a ser el 1, y este apoyo no prescribe el 1
                    % Desprescribimos el 1.
                  else
                    % Queda prescrito todo.
                    % El 1 ya esta
                    % El 2:
                    ES.gdlfij(2*nod) = 1;
                  end
                else % El que estaba prescrito era el 2
                  if theta~=0 && theta~=180 % A menos que los ejes no cambien, solo evenutlamente de signo
                    % Hay que prescribir todo
                    % el 2 ya esta
                   ES.gdlfij(2*nod-1) = 1;
                  end %No aplica else. Si el 2 estaba restringido, va a seguir y el 1 no lo estaba y no debe seguir
                end
             end % No aplica un else pues si ya estaban prescritos... que sigan
           else 
              ES.gdlfij( 2*nod ) = 1; % Se prescribe el segundo grado de libertad solamente 
          end 

          % Modificacion de Fext tradicional 

          ES.Fext( (2*nod) + (-1:0) ) = Tr' * ES.Fext( 2*nod + (-1:0) );

          % Modificacion de Ktotal tradicional

          Kaux = ES.Ktotal ; % Se la define como una nueva variable.
          % Esto es para no andar modificando mucho en la variable de la struct que suele ser mas lento

          % Rotacion de la columna
          for j = 1 : (ES.Nnodo +  ES.NGLX/2) % Como se rotan de a matrices 2x2, el numero de matrices a rotar es la mitad de las filas

            Kaux( (j*2-1):(j*2) , 2*nod + (-1:0)  ) = Kaux( (j*2-1):(j*2) , 2*nod + (-1:0) ) * Tr ;

          end

          % Rotacion de la fila
          for j = 1 : (ES.Nnodo +  ES.NGLX/2) % Como se rotan de a matrices 2x2, el numero de matrices a rotar es la mitad de las filas

            Kaux( 2*nod + (-1:0) , (j*2-1):(j*2) ) = Tr' * Kaux( 2*nod + (-1:0) , (j*2-1):(j*2) );

          end   

          % Actualizacion de K

          ES.Ktotal = Kaux; 

          % Actualizacion de GdLibRot

          ES.GdLiRot( nod ) = ES.CB.Dir.Lin.Desl(i,3);
       end

      if  SeRestringeGX
        nbx = ES.NodAriX( nod ,ArX(AriX))' + (0:1); % Grados de libertad extendidos de este nodo

        % Aca si hay que jugar con  rotaciones

        theta = ES.CB.Dir.Lin.Desl(i,3) - ES.GdLiRot( nbx(2)/2 ); % Angulo de los nuevos ejes respecto al anterior 
                                                        % El e_1 de la nueva base es en el cual desliza.
        if theta==0 || abs(theta)==180 %Caso en el cual desliza segun el eje e_1
           % Modificacion de Upres
           % ---------------------
           ES.Upress( nbx(2) ) = 0;
           % Modificacion de gdlfig
           % ----------------------
           ES.gdlfij( nbx(2) ) = 1;
        elseif abs(theta)==90 || abs(theta)==270
           % Modificacion de Upres
           % ---------------------
           ES.Upress( nbx(1) ) = 0;     
           % Modificacion de gdlfig
           % ----------------------
           ES.gdlfij( nbx(1) ) = 1;
       else
            Tr = [ cosd(theta) , -sind(theta) ;
                 sind(theta) , cosd(theta) ]; %Matriz de rotacion

            Tr = sparse(Tr) ; %Pasada a dispersa para que Matlab/Octave pueda hacer la cuenta de rotacion de la matriz     


             %Modificacion de Upres extendido


             ES.Upres( nbx ) = Tr' * ES.Upres( nbx );

             % Se prescribe el que esta prescripto
             ES.Upres( nbx(2) ) = 0;

             % Modificacion de gdlfij extendido
             % Ver justificacion del if en deslizantes nodales

             if ES.gdlfij( nbx(1) ).^2 + ES.gdlfij( nbx(2) ).^2 ~= 0 % Deben ser los dos nulos (no prescritos) para no entrar en este if
                 if ES.gdlfij( nbx(1) )*ES.gdlfij( nbx(2) ) ==0 % Si alguno de los dos no estaba prescrito
                   if ES.gdlfij( nbx(1) )==1 % El primero estaba prescrito
                     if theta==90 || theta==270 % Si el primero ahora es el segundo (y este es el que queremos prescribir
                      ES.gdlfij( nbx(2) ) = 1;
                      ES.gdlfij( nbx(1) )=0;
                      % Como el 1 original pasa a ser el 2
                      % Lo dejamos prescrito prescribiendo el 2
                      %
                      % Como el 2 original (libre) pasa a ser el 1, y este apoyo no prescribe el 1
                      % Desprescribimos el 1.
                     else
                      % Queda prescrito todo.
                      % El 1 ya esta
                      % El 2:
                      ES.gdlfij( nbx(2) ) = 1;
                     end
                  else % El que estaba prescrito era el 2
                     if theta~=0 && theta~=180 % A menos que los ejes no cambien, solo evenutlamente de signo
                       % Hay que prescribir todo
                       % el 2 ya esta
                       ES.gdlfij( nbx(1) ) = 1;
                     end %No aplica else. Si el 2 estaba restringido, va a seguir y el 1 no lo estaba y no debe seguir
                  end
                end % No aplica un else pues si ya estaban prescritos... que sigan
             else 
                ES.gdlfij( nbx(2)) = 1; % Se prescribe el segundo grado de libertad solamente 
             end 

             % Modificacion de Fext tradicional 

             ES.Fext( nbx ) = Tr' * ES.Fext( nbx );

             % Modificacion de Ktotal tradicional

             Kaux = ES.Ktotal ; % Se la define como una nueva variable.
             % Esto es para no andar modificando mucho en la variable de la struct que suele ser mas lento

             % Rotacion de la columna
             for j = 1 : (ES.Nnodo +  ES.NGLX/2) % Como se rotan de a matrices 2x2, el numero de matrices a rotar es la mitad de las filas

               Kaux( (j*2-1):(j*2) , nbx  ) = Kaux( (j*2-1):(j*2) , nbx ) * Tr ;

             end

             % Rotacion de la fila
             for j = 1 : (ES.Nnodo +  ES.NGLX/2) % Como se rotan de a matrices 2x2, el numero de matrices a rotar es la mitad de las filas

               Kaux( nbx , (j*2-1):(j*2) ) = Tr' * Kaux( nbx , (j*2-1):(j*2) );

             end   

             % Actualizacion de K

             ES.Ktotal = Kaux; 

             % Actualizacion de GdLibRot

             ES.GdLiRot( nbx(2)/2 ) = ES.CB.Dir.Lin.Desl(i,3);        
         end
      end
       
  
   end % Terminar de reocrrer los dos nodos
   
end % Se termina de recorrer la matriz




% ==========================================================================
% === Modificacion de Fext por los prescritos no nulos =====================
% ==========================================================================

ES.Fext = ES.Fext - ES.Ktotal*ES.Upres; % La forma de imponer los prescritos es modificando el vector de fuerza





end