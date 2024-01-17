function ES = Optimizacion(ES,M,c,alpha,TolDeltaV,DeltaVMax,TolM,MaxSit,SmoothingCoef)
% Rutina para intentar optimizar la estructura con el metodo de penalidad
% para las tensiones y lagrangiano aumentado para el area final
% 
% Devuelve ES.psi la estructura optima
%
% Parametros de entrada:
%
% ES la estructura con todas las CB's
% M: Fraccion de volumen buscada
% c: Valor para el factor de penalidad de la diferencia del volumen al
% cuadrado
% alpha: Coeficiente de penalidad para el funcional de tensiones
% TolDeltaV: Minima Fraccion de cambio de volumen con Kappa=1 tal que se
%            considere estructura optima para el Lagrangiano actual
% DeltaVMax: Cuanto es el DeltaVolumen maximo admisible entre iteraciones
% TolM: Fraccion de volumen respecto a la cual M se considera que se llego
% MaxSit: Maximo de subiteraciones en la busqueda lineal. Si se superan se
%         modifica el lagrangiano aunque no sea tan correcto

% ==============================================================================
% === Creacion de una carpeta donde guardar figuras ============================
% ==============================================================================

CurrTime=clock;
FigFolder=['Figuras_Prueba_Del_',num2str(CurrTime(3)),'-',num2str(CurrTime(2)),'-',num2str(CurrTime(1)),'_',num2str(CurrTime(4),'%02.f'),num2str(CurrTime(5),'%02.f'),num2str(floor(CurrTime(6)),'%02.f')];

mkdir(FigFolder)

% Creacion del informe en texto:

%FID=fopen([FigFolder,'/ReporteDeIteraciones.txt'],'w');


% ==============================================================================
% === Matriz de Masa ===========================================================
% ==============================================================================

Masa = Masa_Mat(ES);
MasaNoMust = Masa(~ES.Must,~ES.Must);

% ==============================================================================
% === Volumen inicial ==========================================================
% ==============================================================================

% Normalizar psi:

normpsi = sqrt( ES.psi'*Masa*ES.psi );

ES.psi = ES.psi/normpsi;
if sum(ES.Must)>0
ES.psi(ES.Must)=abs(ES.psi(ES.Must)); % Desde ya los obligatorios los hacemos presentes
end
% Calculo del volumen completo

ESAux=ES; ESAux.psi=ones(length(ES.psi),1);
ESAux=interfase(ESAux); % Chequeo inicial de interfases. Es posible que haya habido cambios y no estan aun actualiados
Omega0 = VolumenXFEM(ESAux);
clear EsAux;

% Volumen inicial

ES=interfase(ES);
Omega=VolumenXFEM(ES);

close all
figure
figXFEM(ES,'Mat',0,0)
% Cositas para el grafico
set(gcf,'Position',[585 226 924 465])
set(gca,'XColor', 'none','YColor','none')
set(gca,'color','none')
xlim([min(ES.Mnodo(:,2)) max(ES.Mnodo(:,2))])
ylim([min(ES.Mnodo(:,3)) max(ES.Mnodo(:,3))])
set(gca,'LooseInset',max(get(gca,'TightInset'),0.02))
title(['Iteracion 0: Volumen ',num2str(100*Omega/Omega0),'%'])
set(gcf, 'InvertHardCopy', 'off'); % setting 'grid color reset' off
axis equal
FigName = [FigFolder,'/Final_Iter0.png'];
set(gcf,'Color',[1 1 1]);
print(FigName,'-dpng')


%fprintf(FID,'\r\n%s\r\n\r\n','Iteracion 0:');
%fprintf(FID,'%s \r\n',['Volumen ',num2str(100*Omega/Omega0),'%']);

% ==============================================================================
% === Iteracion de optimizacion ================================================
% ==============================================================================

n = 1; % Numero de iteracion actual

Lambda = 0; %Primer valor para el Lagrangiano del termino lineal del volumen

flagV = 1; %Una flag para cortar no solo por Tolerancia M sino por optimalidad

while flagV

   DTAcum = zeros(ES.Nnodo,1); %Variable que guarda la Derivada topologica total 
   Fomega = 0; % Idem para el costo

   for LC=1:ES.NLC % Todos los casos de carga
   % SE COMIENZA CALCULANDO LA ESTRUCTURA ACTUAL

       % Se setean las fuerzas de este caso
       if isfield(ES.CB,'NeuCell') % Se chequea si estan en formato cell. Si no lo estan solo puede haber un load case
           ES.CB.Neu.Punt=ES.CB.NeuCell.Punt{LC};
           ES.CB.Neu.Bord=ES.CB.NeuCell.Bord{LC};
           ES.CB.Neu.Vol=ES.CB.NeuCell.Vol{LC};
       end
       ES = xfem(ES);

       % FUNCION DE COSTO ACTUAL - Complacencia

       Fomega = Fomega + ES.U'* ES.Ktotal * ES.U  ;
       
       % DERIVADA TOPOLOGICA - Complacencia

       ES = DerTop_Complacencia_XFEM(ES);

       DTAcum = DTAcum + ES.DTC;
       
       if alpha>0 % Si hay que chequear tensiones
           % FUNCION DE COSTO ACTUAL - Tensiones

           FuncTension = Funcional_Tension(ES); %Calculo del funcional

           Fomega = Fomega + alpha * FuncTension; % Aumento a la funcion de costo

           % DERIVADA TOPOLOGICA - Tensiones

           ES = DerTop_Tensiones_XFEM(ES); % Calculo de la derivada

           DTAcum = DTAcum + alpha * ES.DTT; % Aumento a la derivada topologica
       end
       
   end

   % Agregado de los lagrangianos
   Fomega = Fomega  + Lambda*(Omega-Omega0*M) + 0.5*c*(Omega-Omega0*M)^2;
   DTMOD = DTAcum - sign(ES.psi)*( Lambda + c * (Omega-Omega0*M) );
   DTMOD = sign(DTMOD).*log(1+abs(DTMOD)); %Suavizado que hace novotny...
   % DIRECCION PARA EL NUEVO LEVEL SET

   Phi = DTMOD.*sign(ES.psi);
   PhiNoMust = Phi(~ES.Must); % Solo nos quedamos con aquellos que no son obligatorios
   PsiNoMust = ES.psi(~ES.Must); %Idem para psi
   
   % ANGULO ENTRE EL LEVEL SET ACTUAL Y ESTA DIRECCION

   normphi = sqrt( PhiNoMust' * MasaNoMust *PhiNoMust );

   theta = acosd ( (PhiNoMust'*MasaNoMust*PsiNoMust) / normphi );

   % SE COMIENZA LA SUBITERACION
   Fold = Fomega; %Valor del funcional que se busca minimizar
   psiOld = ES.psi(~ES.Must); % Valor original de la funcion level set
   Fnew = Fold+1; % Valor del funcional nuevo. Lo comenzamos arbitrariamente
                  % Mayor al anterior y la iteracion termina cuando este
                  % decrece.
   Kappa = 1; % Factor de busqueda lineal
   OmegaOld = Omega/Omega0; % Volumen relativo anterior, para los chequeos de
                            % volumen para estructura optima y restirccion
                            % del Delta Vol max.

   % Nuevo valor para el level set. (Es exactamente la funcion Phi normalizada)
   ES.psi(~ES.Must) = 1/sind(theta) * ( sind( (1-Kappa) * theta ) *psiOld + sind ( Kappa*theta ) * PhiNoMust/normphi  );
   % De esta forma los level set de los que son obligatorios no cambia (al
   % inicio se hicieron positivos)
   ES=interfase(ES); % Nuevo chequeo de interfases
   Omega = VolumenXFEM(ES);  % Nuevo volumen


   DeltaVolMax=abs(Omega/Omega0-OmegaOld); % Fraccion de cambio de volumen maxima
   % en esta iteracion.
   
   TolDeltaVMod = min(TolDeltaV, abs(OmegaOld-M)); % Logicamente si la tolerancia para estructura "optima"
                                                       % Es mayor a la lo que falta para el objetivo.
                                                       % Usamos como tolerancia lo que falta
   
  
  
   if abs(OmegaOld-M)<=TolM && DeltaVolMax<TolDeltaV % Es optima pero ademas se llego al volumen
        flagV=0; % Es el final
        ES.psi(~ES.Must) = psiOld;
         ES=interfase(ES); % Nuevo chequeo de interfases
       Omega = VolumenXFEM(ES);  % Nuevo volumen
   else
                                                       
        if DeltaVolMax<TolDeltaVMod % La estructura ya es muy optima pero no cumple volumen


           ES.psi(~ES.Must)=psiOld;
           ES=interfase(ES);
           Omega=VolumenXFEM(ES);
           while DeltaVolMax<TolDeltaVMod

               Lambda = Lambda + c * (Omega-Omega0*M); % Nuevo lagrangiano
               c=c*2;

               DTMOD = DTAcum - sign(ES.psi)*( Lambda + c * (Omega-Omega0*M) );
               Phi = DTMOD.*sign(ES.psi);
               PhiNoMust = Phi(~ES.Must); % Solo nos quedamos con aquellos que no son obligatorios
               PsiNoMust = ES.psi(~ES.Must); %Idem para psi

               % ANGULO ENTRE EL LEVEL SET ACTUAL Y ESTA DIRECCION

               normphi = sqrt( PhiNoMust' * MasaNoMust *PhiNoMust );

               theta = acosd ( (PhiNoMust'*MasaNoMust*PsiNoMust) / normphi );
               ES.psi(~ES.Must) = 1/sind(theta) * ( sind( (1-Kappa) * theta ) *psiOld + sind ( Kappa*theta ) * PhiNoMust/normphi  );

               ES=interfase(ES); % Nuevo chequeo de interfases
               Omega = VolumenXFEM(ES);  % Nuevo volumen

               DeltaVolMax=abs(Omega/Omega0-OmegaOld); % Fraccion de cambio de volumen maxima
               % en esta iteracion.

               ES.psi(~ES.Must) = psiOld; %Volvemos a la anterior
               ES=interfase(ES);
               Omega=VolumenXFEM(ES);
           end

       else % Este lagrangiano aun puede mejorar la estructura

           DeltaVol = abs(Omega/Omega0-OmegaOld); % Inicializamos el parametro de cambio de volumen

           while DeltaVol > DeltaVMax %Hasta que se alcance la restriccion de volumen
               % se hace sin XFEM para hacerla mas rapido
               Kappa = Kappa/2;
               ES.psi(~ES.Must) = 1/sind(theta) * ( sind( (1-Kappa) * theta ) *psiOld + sind ( Kappa*theta ) * PhiNoMust/normphi  );

               ES=interfase(ES); % Chequeo de interfases
               Omega = VolumenXFEM(ES);

               DeltaVol=abs(Omega/Omega0-OmegaOld);
           end

           % Alcanzada la restriccion de cambio de volumen se busca mejorar el
           % funcional
           sit=0;
           while Fnew > Fold && sit<MaxSit
               sit = sit+1;
               ES.psi(~ES.Must) = 1/sind(theta) * ( sind( (1-Kappa) * theta ) *psiOld + sind ( Kappa*theta ) * PhiNoMust/normphi   );


               ES=interfase(ES); % Chequeo inicial de interfases. Es posible que haya habido cambios y no estan aun actualiados
               Omega = VolumenXFEM(ES);

               Fnew = 0; % Hay que obtenerlo
               Comp = 0; % Complacencia para el reporte
               for LC=1:ES.NLC % Todos los casos de carga
               % SE COMIENZA CALCULANDO LA ESTRUCTURA ACTUAL

                   % Se setean las fuerzas de este caso
                   if isfield(ES.CB,'NeuCell') % Se chequea si estan en formato cell. Si no lo estan solo puede haber un load case
                       ES.CB.Neu.Punt=ES.CB.NeuCell.Punt{LC};
                       ES.CB.Neu.Bord=ES.CB.NeuCell.Bord{LC};
                       ES.CB.Neu.Vol=ES.CB.NeuCell.Vol{LC};
                   end
                   ES = xfem(ES);

                   % FUNCION DE COSTO - Complacencia

                   Fnew = Fnew + ES.U'* ES.Ktotal * ES.U  ; 
                   Comp = Comp + ES.U'* ES.Ktotal * ES.U ;
                   if alpha>0
                       % FUNCION DE COSTO ACTUAL - Tensiones

                       FuncTension = Funcional_Tension(ES); %Calculo del funcional

                       Fnew = Fnew + alpha * FuncTension; % Aumento a la funcion de costo
                   end

               end

               % Agregado de los lagrangianos
               Fnew = Fnew  + Lambda*(Omega-Omega0*M) + 0.5*c*(Omega-Omega0*M)^2;

               % Figura:

               clf
               figXFEM(ES,'Mat',0,0)
               % Cositas para el grafico
               set(gcf,'Position',[585 226 924 465])
               set(gca,'XColor', 'none','YColor','none')
               set(gca,'color','none')
               xlim([min(ES.Mnodo(:,2)) max(ES.Mnodo(:,2))])
               ylim([min(ES.Mnodo(:,3)) max(ES.Mnodo(:,3))])
               axis equal
               set(gca,'LooseInset',max(get(gca,'TightInset'),0.02))               
               set(gcf,'Color',[1 1 1]);
               title(['Iteracion ',num2str(n), '; Subiteracion ',num2str(sit),'; Volumen ',num2str(100*Omega/Omega0),'; Lambda: ',num2str(Lambda),'% c: ',num2str(c),'; \Delta Vmax: ',num2str(100*DeltaVolMax),'%; F_{old}: ',num2str(Fold),'; F_{new}: ',num2str(Fnew)])
               pause(0)


               Kappa = Kappa/2;

           end

           if sit==MaxSit && Fnew>Fold
                Lambda = Lambda + c * (Omega-Omega0*M); % Nuevo lagrangiano
                c=c*2;
           end

           if SmoothingCoef>0
               % Paso de suavizar
               if ES.MeshReg && ES.Nnodo == (ES.Nx+1)*(ES.Ny+1) %Caso sencillo y rapido, todo el rectangulo tiene elementos y la malla es estructurada
                   PSI = reshape(ES.psi,ES.Nx+1,ES.Ny+1); % Forma matricial para el level set
                   [PSI,~]=smoothn(PSI,SmoothingCoef);
                   ES.psi=reshape(PSI,ES.Nnodo,1);
               elseif ES.MeshReg % Al menos es estructurada
                   xPSI = ES.xm:ES.dx:ES.xM;
                   yPSI = ES.ym:ES.dy:ES.yM;
                   [xPSI,yPSI] = meshgrid(xPSI,yPSI);
                   F = scatteredInterpolant(ES.Mnodo(:,2),ES.Mnodo(:,3),ES.psi,'linear');
                   PSI = F(xPSI,yPSI);
                   [PSI,~]=smoothn(PSI,SmoothingCoef);
                   F = scatteredInterpolant(xPSI(:),yPSI(:),PSI(:));
                   ES.psi = F(ES.Mnodo(:,2),ES.Mnodo(:,3));
               else % Este es el caso mas lento de smoothing
                   f = fit( ES.Mnodo(:,2:3), ES.psi,'lowess','Span',SmoothingCoef);
                   ES.psi = f(ES.Mnodo(:,2:3));

               end
           end


           if sum(ES.Must)>0 % Si hay algun nodo obligatorio
               ES.psi(ES.Must) = abs ( ES.psi(ES.Must) ); % a pesar dle smoothing obligamos que los must se mantengan
           end

           % Corregimos el volumen
           ES=interfase(ES); % Chequeo inicial de interfases. Es posible que haya habido cambios y no estan aun actualiados
           Omega = VolumenXFEM(ES);

            % Figura:
            clf
           figXFEM(ES,'Mat',0,0)
           % Cositas para el grafico
           set(gcf,'Position',[585 226 924 465])
           set(gca,'XColor', 'none','YColor','none')
           set(gca,'color','none')
           xlim([min(ES.Mnodo(:,2)) max(ES.Mnodo(:,2))])
           ylim([min(ES.Mnodo(:,3)) max(ES.Mnodo(:,3))])
           set(gca,'LooseInset',max(get(gca,'TightInset'),0.02))
           title(['Iteracion ',num2str(n), '; Volumen ',num2str(100*Omega/Omega0),'% \lambda: ',num2str(Lambda),', c: ',num2str(c)])
           FigName = [FigFolder,'/Final_Iter',num2str(n), '.png'];
           axis equal
           set(gcf, 'InvertHardCopy', 'off'); % setting 'grid color reset' off
           set(gcf,'Color',[1 1 1]);
           print(FigName,'-dpng')


           if rem(n,1)==0 %El numero de iteraciones coincide con este valor
                save([FigFolder,'/Estructura_Iter',num2str(n),'.mat'],'ES','-v7.3')
           end


        %    fprintf(FID,'\r\n%s\r\n\r\n',['Iteracion ',num2str(n),':']);
      %      fprintf(FID,'%s \r\n',['Volumen ',num2str(100*Omega/Omega0),'%']);
      %      fprintf(FID,'%s \r\n',['Lambda ',num2str(Lambda)]);
      %      fprintf(FID,'%s \r\n',['Theta ',num2str(theta),'�']);
       %     fprintf(FID,'%s \r\n',['Complacencia ',num2str(Comp)]);


            n=n+1;
        end
   end

    if abs(Omega/Omega0-M)<=TolM % Una vez que se llegue al volumen se disminuye el suavizado
        SmoothingCoef = SmoothingCoef*1;
    end

end

