function DerTop = AuxDTT(Eps,Sigma,Ea,la,mu,alpha1,alpha2,gamma,CharFunDom,SigmaM,CharFunEst,signoLevelSet,FVol,desplaAdj,InteInsolMat)
% Funcion que a partir de las entradas que se ennumeran a continuacion 
% devuelve el valor de la derivada topologica del funcional de tensiones
% segun Amstutz 2010 con las modificaciones de Delgado 2024 (Tesis)
% 
% Entradas:
% Eps: Vector de deformaciones del elemento en el punto que se este
% considerando. Es un vector de 3 entradas con 
% [\varepsilon_{xx}; \varepsilon_{yy}; \gamma_{xy}];
%
% Sigma: Vector de tensiones del elemento en el punto que se este
% considerando. Es un vector de 3 entradas con 
% [\sigma_{xx}; \sigma_{yy}; \tau_{xy}];
%
% Ea: Con formato igual a Eps, deformaciones del problema adjunto.
%
% la, mu: Parametros de lame. - Se podrian
% calcular pero... para que?
% 
% alpha1 y alpha2, los parametros materiales combinads
%
% gamma: Contraste
%
% CharFunDom: Funcion caracteristica del dominio de tension
%
% SigmaM: Tension maxima.
% 
% CharFunEst: Funcion caracteristica de si el punto es estructura o no
%
% signoLevelSet: Signo de la funcion de nivel
%
% FVol: Densidad de fuerzas de volumen
%
% desplaAdj: "Desplazamiento" del adjunto
% 
% Devuelve:
%
% DerTop: La derivada topologica

% =========================================================================
% === Inicio de la Derivada ===============================================
% =========================================================================

DerTop = 0; % Se inicializan como nulos

% =========================================================================
% === Definicion de tensores ==============================================
% =========================================================================

ICuatro = eye(3); % Tensor identidad de cuarto orden
IDos = [1,1,0;1,1,0;0,0,0]; % Identidad segundo orden tensorial identidad segundo orden

BTensor = 6*mu*ICuatro + (la-2*mu)*IDos;

BTildeTensor = 3*ICuatro - IDos;

TTensor = ((1-gamma)/(1+alpha2*gamma)) * (alpha2*ICuatro + 0.5*((alpha1-alpha2)/(1+alpha1*gamma))*IDos);

STensor = ( (pi * (1-gamma)^2) / (4*SigmaM^2*(1+alpha2*gamma)^2)) * ( 10*ICuatro + ( 3 * ( (1+alpha2*gamma)/(1+alpha1*gamma))^2 - 5)*IDos); 

% =========================================================================
% === Primer termino ======================================================
% =========================================================================

ITS = (ICuatro+TTensor)*Sigma; 
DerTop = DerTop + (gamma-1) * ITS' * Ea; % El producto escalar funciona asi porque la tercer entrada de Ea es el doble que la componente real del tensor.

% =========================================================================
% === Segundo termino =====================================================
% =========================================================================

DerTop = DerTop - (gamma-1) * FVol*desplaAdj; 

% =========================================================================
% === Tercer termino ======================================================
% =========================================================================


SVM2 = ( (Sigma(1) + Sigma(2))^2 - 3 *(Sigma(1)*Sigma(2) - Sigma(3)^2)); % Tension de von mises al cuadrado

t = SVM2/SigmaM^2; % Tension relativa al cuadrado

% Si t es tamaño "razonable" no hay problemas numéricos, pero si t pasa a
% ser muy grande la precisión es baja por lo cual se estudian los
% siguientes casos:
if t<= 1e7 % Se experimento y parecía funcionar hasta t = 1e9 pero para estar seguros
    dPhidt =  (1 + t^32) ^ (-31/32) * t^31;
    PhiPosta = (1+t^32)^(1/32)-1;
else % si se pasa, ya conviene asumir 1+t^32 = t^32 en cuyo caso
    dPhidt = 1;
    PhiPosta = t-1;
end % Recordatorio esto es la derivada de la funcion psi, no de ella respecto a la tension

TBS = TTensor * BTensor  * Sigma;

DerTop = DerTop -CharFunDom * signoLevelSet * dPhidt * TBS'*Eps / SigmaM^2 ;

% =========================================================================
% === Cuarto Termino ======================================================
% =========================================================================

DerTop = DerTop - CharFunDom*signoLevelSet * PhiPosta;

% =========================================================================
% === Quinto Termino ======================================================
% =========================================================================

BTildeS = BTildeTensor*Sigma;
TBTildeS = TTensor * BTildeS;
TBSS = TBTildeS(1)*Sigma(1) + TBTildeS(2)*Sigma(2) + 2*TBTildeS(3)*Sigma(3);
t2 = TBSS/SigmaM^2;
TBTildeTS = TTensor * BTildeTensor * TTensor * Sigma;
TBTSS = TBTildeTS(1)*Sigma(1) + TBTildeTS(2)*Sigma(2) + 2*TBTildeTS(3)*Sigma(3);
t3 = TBTSS/(2*SigmaM^2);
t4 = t + t2 + t3;

if t4<= 1e7 % Se experimento y parecía funcionar hasta t = 1e9 pero para estar seguros
    PhiComb = (1+t4^32)^(1/32)-1;
else % si se pasa, ya conviene asumir 1+t^32 = t^32 en cuyo caso
    PhiComb = t4-1;
end % 


DerTop = DerTop+ CharFunDom * (1-CharFunEst) * (PhiComb - PhiPosta - t2*dPhidt);

% =========================================================================
% === Sexto Termino =======================================================
% =========================================================================

SumaTen= Sigma(1)+Sigma(2); % Suma de las tensiones principales
detS = Sigma(1)*Sigma(2)-Sigma(3)^2;
DifTen = sqrt( SumaTen^2 - 4*detS); % Diferencia de las tensiones principales
SumaTen = SumaTen/SigmaM; % Normalizamos
DifTen = DifTen/SigmaM; % Normalizamos

if SumaTen<-5
    izq = 1;
    pesoi = 0;
elseif SumaTen>=5
    izq = 1000; 
    pesoi = 1;
else 
    izq = floor(1 + (SumaTen + 5)/0.01);
    pesoi = min( (SumaTen - 0.01*floor(SumaTen/0.01))/0.01,1); 
end

if DifTen>=5
    aba = 500;
    pesoaba = 1;
else
    aba = floor( 1 + (DifTen)/0.01); 
    pesoaba = min ( (DifTen - 0.01*floor(DifTen/0.01))/0.01,1);
end

XiS = (1-pesoi) * ( (1-pesoaba)*InteInsolMat(izq,aba) + pesoaba*InteInsolMat(izq,aba+1)) + ...
    pesoi * ( (1-pesoaba)*InteInsolMat(izq+1,aba) + pesoaba*InteInsolMat(izq+1,aba+1));
    
DerTop = DerTop + CharFunDom * CharFunEst * XiS / pi;
% =========================================================================
% === Septimo termino =====================================================
% =========================================================================

SSigma = STensor*Sigma;
SS = SSigma(1) * Sigma(1) + SSigma(2)*Sigma(2) + 2* SSigma(3)*Sigma(3);

DerTop = DerTop + CharFunDom*CharFunEst*dPhidt*SS/pi;