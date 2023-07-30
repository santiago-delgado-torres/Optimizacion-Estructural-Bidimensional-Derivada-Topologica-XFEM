function ES = xfem(ES)
%
% Esta funcion es quien resuelve el problema de 
% elasticidad plana con el metodo de los elementos finitos
% extendido dada una estructura indicada por la variable ES
% que debe ser definida adecuadamente como se muestra en el problema
% de ejemplo. 
%
% Cuando la funcion level set (PSI) es mayor a 0 se considera que es
% se tiene material estrucutral (Material 2 segun ES.PropMat) y si es menor
% a 0 se tiene material vacio (Material 1 segun ES.PropMat).
%

% =========================================================================
% === Deteccion de interfases =============================================
% =========================================================================

ES=interfase(ES);

% =========================================================================
% === Matriz de rigidez ===================================================
% =========================================================================

ES = matrizXFEM(ES); 

% =========================================================================
% === Armado de Vector de Fuerzas =========================================
% =========================================================================

ES = NeumannXFEM(ES);

% =========================================================================
% === Agregado apoyos elasticos ===========================================
% =========================================================================

ES = RobinXFEM(ES); 

% =========================================================================
% === Prescripcion de Desplazamientos =====================================
% =========================================================================

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

ES.U = zeros(2*ES.Nnodo + ES.NGLX , 1); %Se inicializa como nulo el vector de grados de libertad solucion
ES.U( ES.gdlfij==1 ) = ES.Upres( ES.gdlfij==1 ) ; % A los restrictos se les asigna su valor prescripto
ES.U(NoRest) = Ured; % y al resto el solucion del sistema de ecuaciones anterior

% Rotar aquellos que haya que rotar

for j = 1 : (ES.Nnodo + ES.NGLX/2) % Como se rotan de a pares, el numero de puntos a rotar  es la mitad de las filas
  
   theta = ES.GdLiRot(j); % Angulo que el nodo esta rotado
   if theta ~= 0
   Tr = [ cosd(theta) , -sind(theta) ;
             sind(theta) , cosd(theta) ]; %Matriz de rotacion

   ES.U( [j*2-1,j*2] ) = Tr*ES.U( [j*2-1,j*2] ); %Se devuelve a la base original
   
   ES.Fext( [j*2-1,j*2] ) = Tr*ES.U( [j*2-1,j*2] ); % Idem fuerzas
   
   % Modificacion de Ktotal tradicional
          
   Kaux = ES.Ktotal ; % Se la define como una nueva variable.
   % Esto es para no andar modificando mucho en la variable de la struct que suele ser mas lento
   
   % Rotacion de la columna
   for k = 1 : (ES.Nnodo + ES.NGLX/2) % Como se rotan de a matrices 2x2, el numero de matrices a rotar es la mitad de las filas
     
     Kaux( (k*2-1):(k*2) , 2*j + (-1:0)  ) = Kaux( (k*2-1):(k*2) , 2*j + (-1:0) ) * Tr' ;
      
   end
   
   % Rotacion de la fila
   for k = 1 : (ES.Nnodo + ES.NGLX/2) % Como se rotan de a matrices 2x2, el numero de matrices a rotar es la mitad de las filas
     
     Kaux( 2*j + (-1:0) , (k*2-1):(k*2) ) = Tr * Kaux( 2*j + (-1:0) , (k*2-1):(k*2) );
     
   end   
        
   % Actualizacion de K
   
   ES.Ktotal = Kaux;    
   end
end   




end