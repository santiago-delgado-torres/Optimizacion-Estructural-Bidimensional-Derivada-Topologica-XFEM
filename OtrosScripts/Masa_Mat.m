function Masa = Masa_Mat(ES)
% Como la matriz de masa va a ser para productos internos entre variables que
% viven en los nodos y no tienen componentes extendidas. No es mas que hacer la integral
% de las funciones de forma tradicionales.
%
% Igual como puede ser una matriz pesada se la hace dispersa.


KPos = 0; %Un indice de posicion que es comodo para ensamblar los vectores
TAM = 3 * 3 * ES.Nelem; % Cada elemento tiene 3 filas y 3 columnas en esta matriz

KVec = zeros(TAM,1); % Valores de las entradas
KLin = zeros(TAM,1); % Numeros de fila
KCol = zeros(TAM,1); % Numeros de columna


for ele = 1:ES.Nelem

    Ke = zeros(3); %Se inicializa
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


    % 3 PUNTOS DE GAUSS EN TRIANGULOS.
    % (3 puntos alcanza para integrar exactamente polinomios de
    % segundo grado como es el caso)

    etaG = [0.5 0 0.5]; % Coordenadas eta de los puntos.
    xiG = [0.5 0.5 0]; % Coordenadas xi de los puntos
    PesoG = [1/6 1/6 1/6]; % Peso de los puntos.


    for puntoG = 1:3

      N1 = 1-etaG(puntoG)-xiG(puntoG);
      N2 = etaG(puntoG);
      N3 = xiG(puntoG);

      N=[N1, N2,N3];

      Ke = Ke + Jo * PesoG(puntoG) * (N'*N);

    end

    DOF=ne'; %Inicializo

    DOFMAT=repmat(DOF,1,3); %Se arma una matriz auxiliar de los grados de libertad para definir la matriz dispersa

    KLin( (KPos+1):(KPos+9) )=DOFMAT(:); % MATRIZ(:) recorre por columnas, entonces efectivamente esto va a devolver las filas donde se debe ubicar KVec.
    DOFMAT=transpose(DOFMAT); % Truco para hallar las columnas

    KCol( (KPos+1):(KPos+9) )=DOFMAT(:);

    KVec( (KPos+1):(KPos+9) )=Ke(:);

    KPos=KPos+9;
end %Fin de la iteracion en los elementos.

Masa = sparse(KLin, KCol, KVec, ES.Nnodo ,ES.Nnodo);