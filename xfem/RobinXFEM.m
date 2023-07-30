function ES = RobinXFEM(ES)

% Rutina para incorporar los efectos de apoyos elasticos nodales.
% Los mismos se modelan modificando la matriz de rigidez
%
% Modifica
%
% Ktotal: Matriz del sistema

% =========================================================================
% === Nodos biapoyados ====================================================
% =========================================================================

for i = 1:size(ES.CB.Rob.Fijo,1) % Recorremos todas las entradas de la matriz asociada
   
   % Estos no requieren ningun tratamiento especial mas que aumentar la rigidez de 2 entradas de la matrizXFEM
   
   % Rigidez segun X
   ES.Ktotal( 2*ES.CB.Rob.Fijo(i,1)-1 , 2*ES.CB.Rob.Fijo(i,1)-1 ) = ES.Ktotal( 2*ES.CB.Rob.Fijo(i,1)-1 , 2*ES.CB.Rob.Fijo(i,1)-1 ) + ES.CB.Rob.Fijo(i,2);
   
   % Rigidez segun Y
   ES.Ktotal( 2*ES.CB.Rob.Fijo(i,1) , 2*ES.CB.Rob.Fijo(i,1)   ) = ES.Ktotal( 2*ES.CB.Rob.Fijo(i,1) , 2*ES.CB.Rob.Fijo(i,1) ) + ES.CB.Rob.Fijo(i,3);
   
end

% =========================================================================
% === Nodos deslizantes ===================================================
% =========================================================================

for i = 1:size(ES.CB.Rob.Desl,1) % Recorremos todas las entradas de la matriz asociada
   
   
   theta = ES.CB.Rob.Desl(i,2); % Angulo segun desliza
   % Como rotar la matriz puede llevar mucho se decide solo rotar si los angulos no son los ortogonales

   if theta==0 || theta==180
     % Si desliza segun la horizontal, el apoyo elastico es en 2*GL 
     ES.Ktotal( 2*ES.CB.Rob.Desl(i,1), 2*ES.CB.Rob.Desl(i,1) ) = ES.Ktotal( 2*ES.CB.Rob.Desl(i,1), 2*ES.CB.Rob.Desl(i,1) ) +  ES.CB.Rob.Desl(i,3);
   elseif theta==90 || theta==270
     % Si desliza segun la vertical el apoyo elastico  es en 2*GL-1
     ES.Ktotal( 2*ES.CB.Rob.Desl(i,1)-1, 2*ES.CB.Rob.Desl(i,1)-1 ) = ES.Ktotal( 2*ES.CB.Rob.Desl(i,1)-1, 2*ES.CB.Rob.Desl(i,1)-1 ) +  ES.CB.Rob.Desl(i,3);
   else % Caso de apoyo inclinado
     % Este lleva un trabajo levemente mayor.
     % Al tener apoyos inclinados. Lo que suele hacerse
     % es rotar las 2 filas y las 2 columnas de la matriz asociadas al desplazamiento
     % del nodo con apoyo inclinado y luego que el vector de fuerzas externas y el 
     % vector de desplazamientos esten expresados en una base solidaria a la inclinacion del apoyo
     
     % Como el apoyo es elastico en esta funcion, solamente modifica la matriz de rigidez
     % Por este motivo primero se rota, se modifica la rigidez y se desrota
   
     Tr = [ cosd(theta) , -sind(theta) ;
            sind(theta) , cosd(theta) ]; %Matriz de rotacion
          
     Tr = sparse(Tr) ; %Pasada a dispersa para que Matlab/Octave pueda hacer la cuenta de rotacion
          
     Kaux = ES.Ktotal ; % Se la define como una nueva variable.
     % Esto es para no andar modificando mucho en la variable de la struct que suele ser mas lento
   
     % Rotacion de la columna
     for j = 1 : (ES.Nnodo + ES.NGLX/2) % Como se rotan de a matrices 2x2, el numero de matrices a rotar es la mitad de las filas
     
       Kaux( (j*2-1):(j*2) , 2*ES.CB.Rob.Desl(i,1) + (-1:0)  ) = Kaux( (j*2-1):(j*2) , 2*ES.CB.Rob.Desl(i,1) + (-1:0) ) * Tr ;
     
     end
   
     % Rotacion de la fila
     for j = 1 : (ES.Nnodo + ES.NGLX/2) % Como se rotan de a matrices 2x2, el numero de matrices a rotar es la mitad de las filas
     
       Kaux( 2*ES.CB.Rob.Desl(i,1) + (-1:0) , (j*2-1):(j*2) ) = Tr' * Kaux( 2*ES.CB.Rob.Desl(i,1) + (-1:0) , (j*2-1):(j*2) );
     
     end
   
     % Ahora el grado de libertad 2*numero - 1 es la componente deslizante y la 2*numero es la elastica
     % Se debe incorporar solo en la segunda
     Kaux( 2*ES.CB.Rob.Desl(i,1) , 2*ES.CB.Rob.Desl(i,1) ) = Kaux( 2*ES.CB.Rob.Desl(i,1) , 2*ES.CB.Rob.Desl(i,1) ) + ES.CB.Rob.Desl(i,3);

     % Desrotacion de la columna
     for j = 1 : (ES.Nnodo + ES.NGLX/2) % Como se rotan de a matrices 2x2, el numero de matrices a rotar es la mitad de las filas
     
       Kaux( (j*2-1):(j*2) , 2*ES.CB.Rob.Desl(i,1) + (-1:0)  ) = Kaux( (j*2-1):(j*2) , 2*ES.CB.Rob.Desl(i,1) + (-1:0) ) * Tr' ;
     
     end   
   
     % Desrotacion de la fila
     for j = 1 : (ES.Nnodo + ES.NGLX/2) % Como se rotan de a matrices 2x2, el numero de matrices a rotar es la mitad de las filas
     
       Kaux( 2*ES.CB.Rob.Desl(i,1) + (-1:0) , (j*2-1):(j*2) ) = Tr * Kaux( 2*ES.CB.Rob.Desl(i,1) + (-1:0) , (j*2-1):(j*2) );
     
     end   
   
     % Actualizacion de K
   
     ES.Ktotal = Kaux;
   
   end
   
   
   
   
end





end
