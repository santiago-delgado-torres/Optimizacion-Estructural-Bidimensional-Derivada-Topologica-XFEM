function ES = MeshNotReg(ES,FileName)
% Rutina para leer una malla no estructurada de un archivo particular.
% De momento solo implementada para archivos t3s generados a partir del
% software BlueKenue 
% https://nrc.canada.ca/en/research-development/products-services/software-applications/blue-kenuetm-software-tool-hydraulic-modellers
% 
% Devuelve o actualiza:
%
% ES.Nnodo:    Cantidad de nodos
% ES.Mnodo:    Matriz de nodos, primera columna indice, segunda coord x,
%              tercera coord y.
% ES.Nelem:    Cantidad de elementos
% ES.Melem:    Matriz de conectividad, primera columna indice de elemento,
%              segunda columna indice de material (en XFEM esta no sirve
%              pero queda incorporada para codigos de FEM), siguientes 3
%              numeros de nodos del elemento.
% 

% =========================================================================
% === Chequeo de archivo ==================================================
% =========================================================================

if strcmp(FileName( (end-3):end ),'.t3s')
    FileType = 1; % Es archivo tipo BlueKenue
else
    error(['No se reconoce el formato de malla ',FileName( (end-3):end ),'. Chequee que sea el tipo correcto'])
end

% =========================================================================
% === Ejecucion para archivo tipo BlueKenue ===============================
% =========================================================================

if FileType==1
    
    FileMesh=fopen(FileName); 
    
    flag=1; subflag1=0; subflag2=0;
    while flag %Este while es hasta que se encuentren dos "keywords"
      Linea=fgetl(FileMesh);
      if length(Linea)>=10 && strcmp(Linea(1:10),':NodeCount') %La keyword de cantidad de nodos.
          subflag1=1;
          NNodes=str2double(Linea(12:end));
      end
      if length(Linea)>=13 && strcmp(Linea(1:13),':ElementCount') %La keyword de cantidad de elementos.
          NElem=str2double(Linea(15:end));
          subflag2=1; 
      end

      flag=1-subflag1*subflag2; %Este numero solamente se hara 0 cuando ambas flags valgan 1, es decir ambas keywords fueron encontradas.
    end

    for i=1:3
        [~]=fgetl(FileMesh); %Avanzamos 3 filas.
    end

    XNodes=zeros(NNodes,1); % Inicializado vector de coordenadas x de la malla
    YNodes=zeros(NNodes,1); %IInicializado vector de coordenadas y de la malla
    Conect=zeros(NElem,3); %Inicializada Matriz de conectividad

    for i=1:NNodes
        Linea=fgetl(FileMesh);
        Valores=str2num(Linea);
        XNodes(i)=Valores(1); YNodes(i)=Valores(2);
    end
    for i=1:NElem
        Linea=fgetl(FileMesh);
        Valores=str2num(Linea);
        Conect(i,:)=Valores;
    end
    
    % Carga a formato del trabajo
    ES.Nnodo = NNodes;
    ES.Nelem = NElem;
    ES.Mnodo = [ (1:NNodes)' , XNodes, YNodes ];
    ES.Melem = [ (1:NElem)', 2*ones(NElem,1), Conect ];

end

