function ES = interfase(ES)

% Rutina para hallar que elementos tienen interfase, que elementos
% son puramente solidos, que elementos son puramente vacios y cuales son
% los grados de libertad extendidos.
%
% Devuelve o actualiza:
%
% ES.EI:     indica si el elemento tiene interfase (Es decir, es extendido)
% ES.VA:     indica si el elemento esta enteramente en el material vacio
% ES.AriX:   lista de aristas extendidas. Son dos columnas con los dos
%            nodos que corresponden.
% ES.NGLX:  indica numero de grados de libertad extendidos
% ES.EGLX:   indica numeracion de aristas extendidas
%            Es una matriz de 2 columas con entradas nulas para los elementos 
%            no extendidos, pero en los extendidos indica los números de
%            aristas, dentro del listado de aristas extendidas, que
%            corresponde.
% ES.NodAriX: Es una matriz dispersa que indica en las filas el numero de
%             nodo, en columna el numero de arista extendida (sería como un
%             indicador del g) y el vaor presentado es el numero de grado
%             de libertad segun X asociado a esta extension
%
% =========================================================================
% === Tolerancia en psi ===================================================
% =========================================================================
ES.psi((ES.psi >= 0) & (ES.psi <  ES.Tolpsi)) =  ES.Tolpsi;
ES.psi((ES.psi <= 0) & (ES.psi > -ES.Tolpsi)) = -ES.Tolpsi;

% =========================================================================
% === Bucle para detectar elementos con interfase y vacios ================
% =========================================================================

ES.EI = zeros(ES.Nelem,1); %Se inicializa el primer vector. Es una variable logica
ES.VA = zeros(ES.Nelem,1); % Idem
ES.AriX = zeros(3*ES.Nelem,2);  % Se inicia con una longitud mayor a la necesaria
ES.NGLX = 0; % Se inicializa nulo y luego se comienza a sumar
ES.EGLX = zeros(ES.Nelem,2); 
NAriX = 0;
% Elementos para la matriz dispersa
FilasX = zeros(3*ES.Nelem*4,1); %Como cada arista tendra entre 3 y cuatro nodos la inicialzamos como 4 veces mas grande a la inicializacion de las aristas
ColumnasX = FilasX;
ValX = FilasX;
NX = 2*ES.Nnodo+1; % Indice de cuales son los grados que se crean
posMat = 1;


for ele = 1:ES.Nelem 
    
    ne = ES.Melem(ele,3:5); % Nodos del elemento
    psie = ES.psi(ne); % Vector de la funcion de nivel en los nodos del elemento
    
    if ( psie(1) * psie(2) < 0 ) || ( psie(1) * psie(3) < 0 ) % Caso en que haya cambio de signo en la funcion de nivel
        % Por Bolzano hay interfas
        
        ES.EI(ele) = 1 ; %Aclaramos que este es un elemento con interfaz (extendido)
        
        pos=1;
        
        if psie(1)*psie(2) < 0
            % la interfase pasa por la arista entre 1 y 2
            
            nA = ne(1); % Primer nodo de la arista
            nB = ne(2); % Idem segundo nodo
            ind2 = find(ES.AriX(1:NAriX,1)==nA | ES.AriX(1:NAriX,2) == nA);
            ind = find(ES.AriX(ind2,1)==nB | ES.AriX(ind2,2) == nB);
            if isempty(ind) % Si es la primera vez que se encuentra esta arista
                NAriX = NAriX+1;
                ES.AriX(NAriX,:)=[nA,nB]; % Se la agrega a la lista
                ES.NGLX = ES.NGLX + 6;  % Esta implica al menos, grados de libertad en los 3 nodos del elemento
                ES.EGLX(ele,pos) = NAriX;
                
                FilasX(posMat:(posMat+2)) = ne'; % Corresponden los 3 nodos
                ColumnasX(posMat:(posMat+2)) = NAriX; % Todos son los de esta linea
                ValX(posMat) = NX; posMat=posMat+1; NX=NX+2;
                ValX(posMat) = NX; posMat = posMat+1; NX=NX+2;
                ValX(posMat) = NX; posMat= posMat+1; NX=NX+2;
                
            else
                ES.NGLX = ES.NGLX + 2; % Como ya existia, solamente agrega grado al nodo del nuevo elemento
                ES.EGLX(ele,pos) = ind2(ind); %ES la arista ya conocida
                
                FilasX(posMat)= ne( ne~=nA & ne~=nB ); %Solo se agrega el que falta
                ColumnasX(posMat) = ind2(ind); % Efecto similar al ES.EGLX
                ValX(posMat) = NX; posMat=posMat+1; NX=NX+2;
            end
            
            pos = pos+1;
            
            
        end
        
        if psie(1)*psie(3) < 0
            % la interfase pasa por la arista entre 1 y 3
            
            nA = ne(1); % Primer nodo de la arista
            nB = ne(3); % Idem segundo nodo
            ind2 = find(ES.AriX(1:NAriX,1)==nA | ES.AriX(1:NAriX,2) == nA);
            ind = find(ES.AriX(ind2,1)==nB | ES.AriX(ind2,2) == nB);
             if isempty(ind) % Si es la primera vez que se encuentra esta arista
                NAriX = NAriX+1;
                ES.AriX(NAriX,:)=[nA,nB]; % Se la agrega a la lista
                ES.NGLX = ES.NGLX + 6;  % Esta implica al menos, grados de libertad en los 3 nodos del elemento
                ES.EGLX(ele,pos) = NAriX;
                
                FilasX(posMat:(posMat+2)) = ne'; % Corresponden los 3 nodos
                ColumnasX(posMat:(posMat+2)) = NAriX; % Todos son los de esta linea
                ValX(posMat) = NX; posMat=posMat+1; NX=NX+2;
                ValX(posMat) = NX; posMat = posMat+1; NX=NX+2;
                ValX(posMat) = NX; posMat= posMat+1; NX=NX+2;
                
            else
                ES.NGLX = ES.NGLX + 2; % Como ya existia, solamente agrega grado al nodo del nuevo elemento
                ES.EGLX(ele,pos) = ind2(ind); %ES la arista ya conocida
                
                FilasX(posMat)= ne( ne~=nA & ne~=nB ); %Solo se agrega el que falta
                ColumnasX(posMat) = ind2(ind); % Efecto similar al ES.EGLX
                ValX(posMat) = NX; posMat=posMat+1; NX=NX+2;
            end
            
            pos = pos+1;
            
            
        end
            
    
        if psie(2)*psie(3) < 0
            % la interfase pasa por la arista entre 2 y 3
            
            nA = ne(2); % Primer nodo de la arista
            nB = ne(3); % Idem segundo nodo
            ind2 = find(ES.AriX(1:NAriX,1)==nA | ES.AriX(1:NAriX,2) == nA);
            ind = find(ES.AriX(ind2,1)==nB | ES.AriX(ind2,2) == nB);
            if isempty(ind) % Si es la primera vez que se encuentra esta arista
                NAriX = NAriX+1;
                ES.AriX(NAriX,:)=[nA,nB]; % Se la agrega a la lista
                ES.NGLX = ES.NGLX + 6;  % Esta implica al menos, grados de libertad en los 3 nodos del elemento
                ES.EGLX(ele,pos) = NAriX;
                
                FilasX(posMat:(posMat+2)) = ne'; % Corresponden los 3 nodos
                ColumnasX(posMat:(posMat+2)) = NAriX; % Todos son los de esta linea
                ValX(posMat) = NX; posMat=posMat+1; NX=NX+2;
                ValX(posMat) = NX; posMat = posMat+1; NX=NX+2;
                ValX(posMat) = NX; posMat= posMat+1; NX=NX+2;
                
            else
                ES.NGLX = ES.NGLX + 2; % Como ya existia, solamente agrega grado al nodo del nuevo elemento
                ES.EGLX(ele,pos) = ind2(ind); %ES la arista ya conocida
                
                FilasX(posMat)= ne( ne~=nA & ne~=nB ); %Solo se agrega el que falta
                ColumnasX(posMat) = ind2(ind); % Efecto similar al ES.EGLX
                ValX(posMat) = NX; posMat=posMat+1; NX=NX+2;
            end
                        
        end
            
    elseif psie(1) < 0 % Si no hay cambio de signo, entonces con chequear que cualquiera
        % es menor a 0 alcanza para decir que es un vacio.
        
        ES.VA(ele) = 1; 
    end
    
    
end

ES.AriX = ES.AriX(1:NAriX,:); %Reduccion de la previamente iniciada

% Matriz dispersa
FilasX = FilasX(1:(posMat-1));
ColumnasX=ColumnasX(1:(posMat-1));
ValX = ValX(1:(posMat-1));

ES.NodAriX = sparse(FilasX,ColumnasX,ValX,ES.Nnodo,NAriX);

end