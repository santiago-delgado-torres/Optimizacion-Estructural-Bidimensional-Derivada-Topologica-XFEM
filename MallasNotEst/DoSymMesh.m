function ES = DoSymMesh(ES,Eje,Pos)
% Rutina para simetrizar una malla no estructurada.
% 
% Devuelve ES.Mnodo, ES.Nnodo, ES.Melem y ES.Nelem actualizados
%
% Parametros de entrada:
%
% ES la estructura con la malla
% Eje: Eje de simetria. Opciones 'x' o 'y'. Eventualmente se puede
% considerar otro caso pero por ahora lo ignoramos. Probablemente sea mas
% facil rotar la malla simetrizar y desrotarla.
% Pos: Coordenada (en el otro eje) del eje de simetria.

% ==============================================================================
% === Seleccion del Eje ========================================================
% ==============================================================================

if strcmp(Eje,'x')
    IndexEje = 2; %Vamos a comparar valores en el eje y
elseif strcmp(Eje,'y')
    IndexEje = 1; % Idem pero en x
else
    error(['La opcion de Eje: ',Eje,' no está habilitada. Debe seleccionar x o y'])
end

% ==============================================================================
% === Simetrizar nodos =========================================================
% ==============================================================================

Correspondientes = zeros(1,ES.Nnodo); % Nodos correspondientes a los simetrizados

indexNodo = ES.Nnodo + 1; % Numero del primer nodo nuevo

NewNodes = zeros(ES.Nnodo,3); %Se inicializa matriz de nodos a agregar

for n = 1:ES.Nnodo
    xN = ES.Mnodo(n,2);
    yN = ES.Mnodo(n,3);
    
    if IndexEje == 2 % Se estudia segun y
        if yN == Pos % Nodo que no hay que simetrizar
            Correspondientes(n)=n;
        else
            Correspondientes(n) = indexNodo;
            xNew = xN;
            yNew = Pos + (Pos-yN); %Desde el punto de vista teorico, es lo mismo que 2Pos - yN pero es mas claro pensar que es la posicion + que tan abajo esta el yN.
            NewNodes(indexNodo-ES.Nnodo,:)=[indexNodo,xNew,yNew];
            indexNodo = indexNodo+1;
        end
    else 
        if xN == Pos % Nodo que no hay que simetrizar
            Correspondientes(n)=n;
        else
            Correspondientes(n) = indexNodo;
            yNew = yN;
            xNew = Pos + (Pos-xN);
            NewNodes(indexNodo-ES.Nnodo,:)=[indexNodo,xNew,yNew];
            indexNodo = indexNodo+1;
        end
    end
end

% Se descartan los que quedaron vacios:
NewNodes=NewNodes(1:(indexNodo-1-ES.Nnodo),:);


% ==============================================================================
% === Simetrizar elementos  ====================================================
% ==============================================================================

NewElem = 0*ES.Melem;

for ele = 1:ES.Nelem
    NewElem(ele,:) = [ele + ES.Nelem , ES.Melem(ele,2), Correspondientes(ES.Melem(ele,5:-1:3))]; % El 5:-1:3 es para mantener la orientacion delos elementos
end

% ==============================================================================
% === Se devuelve la malla simetrizada =========================================
% ==============================================================================
ES.Nnodo=indexNodo-1;
ES.Mnodo=[ES.Mnodo;NewNodes];
ES.Nelem=2*ES.Nelem;
ES.Melem=[ES.Melem;NewElem];

