function DerTop = AuxDTT(Eps,Sigma,Ea,lambda,mu,rho,eta,xi,gamma1,gamma0,CharFun,SigmaM)
% Funcion que a partir de las entradas que se ennumeran a continuacion 
% devuelve el valor de la derivada topologica del funcional de tensiones
% segun Amstutz 2010.
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
% lambda y mu: Parametros de Lame.
% 
% rho, eta y xi son los valores de rho, eta y xi segun el paper.
% 
% gamma1 y gamma0 son los valores de los coeficientes gamma del paper,
% siendo el primero el del material a agregar y el segundo el del actual.
%
% CharFun: Funcion caracteristica del dominio. Con fines practicos es ver
% si en los nodos que se va a agregar la derivada se integra el funcional o
% no.
%
% SigmaM: Tension maxima.
% 
% Devuelve:
%
% DerTop: La derivada topologica

% =========================================================================
% === Inicio de la Derivada ===============================================
% =========================================================================

DerTop = 0; 

% =========================================================================
% === Parametros que importan a lo largo de toda la derivada ==============
% =========================================================================

BTens = [ 4*mu + lambda ,   lambda - 2*mu  ,  0   ;
      lambda - 2*mu ,   4*mu + lambda  ,  0   ;
            0       ,        0         , 6*mu];

BtildeTens = [2,-1,0;
              -1,2,0;
              0,0,3];

BS = BTens*Sigma; 

SVM2 = 0.5 * ( BS(1)*Eps(1) + BS(2)*Eps(2) + BS(3)*Eps(3));

t = SVM2 / (SigmaM)^2; 
if t<= 1e7 
    dPhidSVM2 =  (1 + t^32) ^ (-31/32) * t^31;
else 
    dPhidSVM2 = 1;
end 

k1 = dPhidSVM2 * CharFun;

gamma = gamma1/gamma0;

T1 = eta * eye(3); 
T2 = 0.5 * ( xi - eta ) / ( 1 + gamma*xi ) ;
T2 = T2 * [ 1 , 1 , 0;
            1, 1, 0;
            0, 0, 0]; 
T = T1+T2; 


% =========================================================================
% === Primer Termino de la Derivada =======================================
% =========================================================================

TBS = T*BS; 
PrimerTermino = - ( gamma1 - gamma0 ) *  ( rho * k1 * ( TBS(1)*Eps(1) + TBS(2)*Eps(2) + TBS(3)*Eps(3) ) );

DerTop = DerTop + PrimerTermino;

% =========================================================================
% === Segundo Termino de la Derivada ======================================
% =========================================================================

rhoTminusI = rho*T - eye(3);
STS = rhoTminusI*Sigma; 
SegundoTermino = - ( gamma1 - gamma0 ) *   ( STS(1)*Ea(1) + STS(2)*Ea(2) + STS(3)*Ea(3) ) ;

DerTop = DerTop + SegundoTermino;

% =========================================================================
% === Tercer Termino de la Derivada =======================================
% =========================================================================

BtildeS = BtildeTens*Sigma;

t1 = 0.5 * (BtildeS(1)*Sigma(1) + BtildeS(2)*Sigma(2) + 2*BtildeS(3)*Sigma(3));

TBtildeS = T*BtildeS; 

t2 = -rho*( TBtildeS(1)*Sigma(1) + TBtildeS(2)*Sigma(2) + 2*TBtildeS(3)*Sigma(3));

TBtildeTS = T*BtildeTens*T*Sigma; 

t3 = -0.5*rho^2 * ( TBtildeTS(1)*Sigma(1) + TBtildeTS(2)*Sigma(2) + 2*TBtildeTS(3)*Sigma(3) );

t4 = (t1+t2+t3)/(SigmaM)^2;


if t4<= 1e7 
    Phi_3 = ( 1 + t4^32 )^(1/32) -1;
else 
    Phi_3 = t4 - 1;
end

TercerTermino = CharFun.*gamma1.*(Phi_3 + rho*k1*( TBtildeS(1)*Sigma(1) + TBtildeS(2)*Sigma(2) + 2*TBtildeS(3)*Sigma(3)));

DerTop = DerTop + TercerTermino;


% =========================================================================
% === Cuarto Termino de la Derivada =======================================
% =========================================================================

trS = Sigma(1)+Sigma(2); 
DetS = Sigma(1)*Sigma(2)-Sigma(3)^2;
SI = 0.5 * ( trS + sqrt( trS^2 -4 * DetS ) ); 
SII = 0.5 * ( trS - sqrt( trS^2 -4*DetS ) ); 

PsiP = PsiPGaussPoints(SI,SII,t1,rho,SigmaM,eta,xi,gamma);

Constante = ( 1 + eta*gamma ) / (1+xi*gamma);

SS = Sigma(1)*Sigma(1) + Sigma(2)*Sigma(2) + 2*Sigma(3)*Sigma(3); 

CuartoTermino = 1/pi * gamma0*CharFun.* ( PsiP + 0.25 * rho^2 * k1 * pi * ( 5*( 2*SS -trS^2 ) +3*Constante^2*trS^2 ) );

DerTop = DerTop + CuartoTermino;

% =========================================================================
% === Quinto Termino de la Derivada =======================================
% =========================================================================

if t<= 1e7
    Phi_5 = ( 1 + t^32 )^(1/32) -1;
else 
    Phi_5 = t - 1;
end

QuintoTermino = - CharFun.*gamma0.*Phi_5;

DerTop = DerTop + QuintoTermino;
