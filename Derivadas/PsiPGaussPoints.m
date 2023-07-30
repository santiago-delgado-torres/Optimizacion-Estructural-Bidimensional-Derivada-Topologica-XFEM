function PsiP = PsiPGaussPoints(SI,SII,SM2,rho,SigmaM,eta,xi,gamma)
% Calculo, mediante integracion en puntos de gauss
% de la funcion psiP del paper de Amstuts 2010.
%
% Entradas:
% SI: Tension principal maxima
% SII: Tension principal minima
% SM2: Tension de VonMises al cuadrado calculada con Btilde
% rho: Constante material determinada previamente
% SigmaM: Tension de von Mises limite
% eta, xi: Los parametros materiales anteriormente definidos
% gamma: El cociente entre gamma1 y gamma0
%
% Devuelve:
% PsiP el valor de ese funcional.

% =========================================================================
% ===== Definicion de los puntos de Gauss =================================
% =========================================================================

etaG = [ 0.932469514203152;
        -0.932469514203152;
         0.661209386466265;
        -0.661209386466265;
         0.238619186083197;
        -0.238619186083197]; % Coordenadas intrinsecas
pesoG = [0.171324492379170;
         0.171324492379170;
         0.360761573048139;
         0.360761573048139;
         0.467913934572691;
         0.467913934572691]; % Pesos de los puntos
     
% =========================================================================
% ===== Integracion en Puntos de Gauss ====================================
% =========================================================================

Jt = 1/2; %Jacobiano para la integral de t
Jtheta = pi/2; %Idem para la integral de theta

PsiP = 0;
Aux = ( 1 + eta*gamma ) / (1+xi*gamma); %Variable auxiliar que aparece muchas veces

for puntoGT = 1:6
    t = 0.5 + 0.5*etaG(puntoGT); % Valor de t en este punto
    ValIntegralTheta = 0; %Valor de la integral en theta
    
    for puntoGA = 1:6
        theta = pi/2 + pi/2*etaG(puntoGA); %Valor de theta en este punto de gauss
        
        % Creacion del Delta
        Delta = rho * (t/2) * ( ( SI^2 - SII^2 )*(2+3*Aux)*cos(theta) + 3*(SI-SII)^2*(2-3*t)*cos(2*theta));
        Delta = Delta + rho^2*(t^2/4) * ( 3 *( SI + SII )^2*Aux^2 + (SI-SII)^2 *( 3 *( 2 - 3*t )^2 + 4*(cos(theta))^2 ) + 6 * Aux * (SI^2-SII^2) * (2-3*t)*cos(theta));
        
        % Primer valor que figura dentro de un phi
        t1 = ( SM2 + Delta ) / (SigmaM^2);
        % Segundo valor que figura dentro de un phi
        t2 = ( SM2 ) / (SigmaM^2);
        
        % Si t es tamaño "razonable" no hay problemas numéricos, pero si t pasa a
        % ser muy grande la precisión es baja por lo cual se estudian los
        % siguientes casos:
        if t2<= 1e7 % Se experimento y parecía funcionar hasta t = 1e9 pero para estar seguros
            Phi2 = ( 1 + t2^32 )^(1/32) -1;
            dPhidSVM2 = (1/SigmaM)^2 * (1 + t2^32) ^ (-31/32) * t2^31;
        else % si se pasa, ya conviene asumir 1+t^32 = t^32 en cuyo caso
            Phi2 = t2-1;
            dPhidSVM2 = (1/SigmaM)^2;
        end
        
        % Idem analisis con t1
        if t1 <= 1e7
            Phi1 = (1+t1^32)^(1/32)-1;
        else
            Phi1 = t1-1;
        end
        
        Integrando = (1/t^2)* ( Phi1 - Phi2 -dPhidSVM2*Delta );
        

        
        ValIntegralTheta = ValIntegralTheta + Jtheta * pesoG(puntoGA) * Integrando;
    end
    
    PsiP = PsiP + ValIntegralTheta * Jt * pesoG(puntoGT); % Agregado ese punto de Gauss
end


    