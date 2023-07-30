
X=ES.Mnodo(:,2);%Coord X de los nodos

Y=ES.Mnodo(:,3); %Coord Y de los nodos

U=ES.U; 

Psi = ES.psi; % Level set de los nodos

%% Plot de desplazamietnos



paso=5e-3;
vecX=-.5:paso:0.4999;
vecY=-.5:paso:0.4999;

xplot=zeros(length(vecX),length(vecY));
yplot=xplot;
UTot=xplot;
VTot=yplot;
UFEM = xplot;
VFEM = xplot;
UXFEM = xplot;
VXFEM = xplot;
SigmaX = xplot;



xm = min(ES.Mnodo(:,2));
ym = min(ES.Mnodo(:,3));
dx = ES.Lx/ES.Nx;
dy = ES.Ly/ES.Ny;

punto=1;

for i = 1:length(vecX)
    for j = 1:length(vecY)
        xplot(punto)=vecX(i);
        yplot(punto)=vecY(j);
        
        indx = 1+floor( (vecX(i)-xm)/dx);
        indy = 1+floor( (vecY(j)-ym)/dy);
        
        ele = 1+2*(indx-1) + (indy-1)*2*ES.Nx;
        
        
        xizq = xm + dx*(indx-1);
        yaba = ym + dy*(indy-1);
        
        if vecY(j)-yaba <= vecX(i) - xizq %Elemento abajo a la derecha
            ind1=indx + (indy-1)*(ES.Nx+1);
            ind2= ind1+1;
            ind3 = ind2 + ES.Nx+1;
            xi = (vecY(j)-yaba)/dy;
            eta = (vecX(i)-xizq)/dx-xi;
        else
            ind1= indx + (indy-1)*(ES.Nx+1);
            ind3=ind1+ES.Nx+1;
            ind2=ind3+1;
            eta=(vecX(i)-xizq)/dx;
            xi=(vecY(j)-yaba)/dy - eta;
            ele = ele+1;
        end
        
        if xi <0
            xi=0;
        end
        if eta<0
            eta=0;
        end
        
        uele = [U(ind1*2-1),U(ind2*2-1),U(ind3*2-1)]';
        vele = [U(ind1*2),U(ind2*2),U(ind3*2)]';
        
        
        N1 = 1-eta-xi;
        N2 = eta;
        N3 = xi;
        N=[N1 N2 N3];
        
        
         ne = ES.Melem(ele,3:5);
        psie = ES.psi(ne); %Funcion de nivel para los nodos del elemento.
            
        psirep = N*psie;
        
        if psirep > 0
            E = ES.Es;
            nu = ES.nus;
        else
            E = ES.Ev;
            nu = ES.nuv;
        end

        if ES.PEL==2
            E = E/(1-nu^2);
            nu = nu / (1-nu);
        end

        C = (E/(1-nu^2))* [ 1 nu 0;
            nu 1 0;
            0 0 0.5*(1-nu)]; % Matriz tensor elastico
        
        % Geometria del elemento:
        Xe = ES.Mnodo(ne,2);
        Ye = ES.Mnodo(ne,3); 

        % Derivadas de las funciones de forma segun coord. intrinsecas
        dN_detachi=[-1 1 0; -1 0 1];

        % Matriz Jacobiana y su determinante
        J = dN_detachi*[Xe,Ye];
        Jo = det(J);

        % Derivadas de las funciones de forma segun x e y.
        dN_dxy = J\dN_detachi; % Se hace con J\ en vez de lo usual de inv(J)* pues 
        % Matlab dice que asi es mas rapido.
        
        Bfem = zeros(3,6);
        Bfem(1,1:2:5)=dN_dxy(1,:);
        Bfem(2,2:2:6)=dN_dxy(2,:);
        Bfem(3,2:2:6)=dN_dxy(1,:);
        Bfem(3,1:2:5)=dN_dxy(2,:);
        
        if ES.EI(ele)
        
            ArX = ES.EGLX(ele,:); %Numero de aristas extendidas de este elemento
            NArX = ES.AriX(ArX,:); %Numero de nodo de esas aristas
            
            % Punto entre nodo 1 y 2. En este xi = 0
            if psie(1)*psie(2)<0 
                etaR1 = psie(1)/(psie(1)-psie(2));

                if ( ne(1)==NArX(1,1) || ne(1)==NArX(1,2) ) && ( ne(2)==NArX(1,1) || ne(2) == NArX(1,2) ) % El punto pasa por la arista 1
                    Ar1R1 = 1; % Indicador de que es el punto donde g1 vale 1
                    Ar2R1 = 0; 
                else
                    Ar2R1 = 1; %Inidcador de que es el punto donde g2 vale 1
                    Ar1R1 = 0;
                end
            else
                etaR1 = 0.5;
                Ar1R1 = 0; Ar2R1 = 0;
            end

            % Punto entre nodo 2 y 3. En este xi = 1-eta
            if psie(2)*psie(3)<0
                etaR2 = psie(3)/(psie(3)-psie(2));
                if ( ne(2)==NArX(1,1) || ne(2)==NArX(1,2) ) && ( ne(3)==NArX(1,1) || ne(3) == NArX(1,2) ) % El punto pasa por la arista 1
                    Ar1R2 = 1; % Indicador de que es el punto donde g1 vale 1
                    Ar2R2 = 0; 
                else
                    Ar2R2 = 1; %Inidcador de que es el punto donde g2 vale 1
                    Ar1R2 = 0;
                end
            else
                etaR2 = 0.5;
                Ar1R2 = 0; Ar2R2 = 0;
            end  

            % Punto entre nodo 3 y 1. En este eta=0.        
            if psie(3)*psie(1)<0
                xiR3 = psie(1)/(psie(1)-psie(3));
                if ( ne(1)==NArX(1,1) || ne(1)==NArX(1,2) ) && ( ne(3)==NArX(1,1) || ne(3) == NArX(1,2) ) % El punto pasa por la arista 1
                    Ar1R3 = 1; % Indicador de que es el punto donde g1 vale 1
                    Ar2R3 = 0; 
                else
                    Ar2R3 = 1; %Inidcador de que es el punto donde g2 vale 1
                    Ar1R3 = 0;
                end
            else
                xiR3 = 0.5;
                Ar1R3 = 0; Ar2R3 = 0;
            end  
            

            % INTEGRACION EN SUBELEMENTOS
            %
            % Obs: Se dice subelementos como abuso de nomeclatura. En realidad
            % es siempre el mismo elemento y son subtriangulos para
            % integracion y definicion de funciones

            etaSubElements = [0      1      0   etaR1;
                             etaR1 etaR2   0   etaR2;
                             0     etaR1 etaR2   0  ];
            % Ordenado en columnas la anterior variable es las coordenadas eta
            % de los subelementos.

            xiSubElements = [0      0       1      0    ;
                             0   1-etaR2   xiR3  1-etaR2;
                           xiR3   0    1-etaR2  xiR3  ];

            % Ordenado en columnas la anterior variable es las coordenadas xi
            % de los subelementos.

            Ar1SE = [0      0      0     Ar1R1;
                     Ar1R1  Ar1R2  Ar1R3 Ar1R2;
                     Ar1R3  Ar1R1  Ar1R2 Ar1R3 ]; % Formato similar a los anteriores para los valores de g1 en los vertices de los subtriangulos

            Ar2SE = [0      0      0     Ar2R1;
                     Ar2R1  Ar2R2  Ar2R3 Ar2R2;
                     Ar2R3  Ar2R1  Ar2R2 Ar2R3 ]; % Formato similar a los anteriores para los valores de g2 en los vertices de los subtriangulos

            sube = zeros(1,4);
            for ko = 1:4
                sube(ko) = inpolygon(eta,xi,etaSubElements(:,ko),xiSubElements(:,ko));
            end
            subEl = find(sube,1);
            
              
            etase = etaSubElements(:,subEl);% Coordenadas eta de este subtriangulo.
            xise = xiSubElements(:,subEl); % Idem para la coordenada xi.
            Ar1 = Ar1SE(:,subEl); % Valores de Ar para g1
            Ar2 = Ar2SE(:,subEl); %Valores de Ar para g2
            
            % Para hallar los valores de eta eval y xieval (los valores de
            % eta y xi en el subtriangulo) para los valores de eta y xi en
            % el triangulo hay que hacer un cambio de base a lo gal 1
            
            vecETA = [etase(2)-etase(1) , xise(2)-xise(1) ]';
            vecXI = [etase(3) - etase(1) , xise(3)-xise(1) ]';
            
            IdSET = [vecETA, vecXI]; %Esta es la matriz que me pasa un vector en la base del subtriangulo a la del triangulo, con origen en etase(1) y xise(1)
            IdTSE = inv(IdSET); % Esta entonces seria la matriz que me pasa una del triangulo con origen etase(1) xise(1) a la del subtriangulo
            
            
            vecTo = [eta - etase(1) , xi - xise(1) ]';
            
            vecAux = IdTSE*vecTo;
            
            etaeval = vecAux(1);
            xieval = vecAux(2);


            N1 = 1-etaeval-xieval; % funcion de interpolacion clasica N1 para el elemento completo en este punto
            N2 = etaeval; %Idem N2
            N3 = xieval; %Idem N3
            Nse = [N1,N2,N3]; % vector de las mismas.
            
            
            % Funciones para definir los g
            g1 = Nse * Ar1;
            g2 = Nse * Ar2;
            
            NX(1:3)= N*g1;
            NX(4:6) = N*g2;
            
            BX=zeros(3,12);
            
            JacSe = dN_detachi*[etase xise]; %Jacobiano del subelemento
            
            Jse = det( JacSe ); % Determinante jacobiano de este subtriangulo en el triangulo intrinseco
            
            dNse_detachi = JacSe\dN_detachi;
            
            dg1_detachi = dNse_detachi*Ar1;
            dg2_detachi = dNse_detachi*Ar2;
            
            dw1_detachi = dN_detachi*g1 + dg1_detachi*N;
            dw2_detachi = dN_detachi*g2 + dg2_detachi*N;
            
            dw1_dxy = J\dw1_detachi;
            dw2_dxy = J\dw2_detachi;
            
            BX(1,1:2:5) = dw1_dxy(1,:);
            BX(2,2:2:6) = dw1_dxy(2,:);
            BX(3,1:2:5) = dw1_dxy(2,:);
            BX(3,2:2:6) = dw1_dxy(1,:);
            
            
            BX(1,7:2:11) = dw2_dxy(1,:);
            BX(2,8:2:12) = dw2_dxy(2,:);
            BX(3,7:2:11) = dw2_dxy(2,:);
            BX(3,8:2:12) = dw2_dxy(1,:);
            
            
            UX = U( reshape( ES.NodAriX(ne,ArX)  , [6,1]) ) ;
            VX = U( reshape( ES.NodAriX(ne,ArX) +1,[6 1]) );
            
            GLXFEM = [UX(1);VX(1);UX(2);VX(2);UX(3);VX(3);UX(4);VX(4);UX(5);VX(5);UX(6);VX(6)];
        else
            NX = zeros(1,6);
            UX = NX';
            VX = NX';
                 
            GLXFEM = 0;
            BX = 0;
        end
        
        GLFEM = [uele(1);vele(1);uele(2);vele(2);uele(3);vele(3)];
        
        
        EpsFEM = Bfem*GLFEM;
        SigmaFEM = C*EpsFEM;
        EpsXFEM = BX * GLXFEM;
        SigmaXFEM = C*EpsXFEM;
       
        UTot(punto)=N*uele + NX*UX;
        VTot(punto)=N*vele+NX*VX;
        
        UFEM(punto) = N*uele;
        VFEM(punto) = N*vele;
        
        UXFEM(punto) = NX*UX;
        VXFEM(punto) = NX*VX;
        
        SigmaX(punto) = SigmaFEM(1) + SigmaXFEM(1);
        
        punto=punto+1;
    end
end

figure
h=surf(xplot,yplot,UTot);
figure
surf(xplot,yplot,VTot);
figure
h=surf(xplot,yplot,UXFEM);

figure
h=surf(xplot,yplot,VXFEM);

figure
surf(xplot,yplot,SigmaX)

STOP

S = 1;
R = 2;
ep = 1;
al = ES.Ev/ES.Es;
g = al;
E = ES.Ev;
v = 0;
muu = E/2/(1+v);
la = v*E/(1-v^2);
close all

AR0 = (R^2*g*la + 2*R^2*g*muu)/(4*(R^2*muu - ep^2*muu + R^2*g*la + R^2*g*muu + ep^2*g*muu));
AR1 = -(R^4*ep^2*g^2*la^2 + 3*R^4*ep^2*g^2*la*muu + 2*R^4*ep^2*g^2*muu^2 - R^4*ep^2*g*la^2 - 3*R^4*ep^2*g*la*muu - 2*R^4*ep^2*g*muu^2 - R^2*ep^4*g^2*la^2 - 3*R^2*ep^4*g^2*la*muu - 2*R^2*ep^4*g^2*muu^2 + R^2*ep^4*g*la^2 + 3*R^2*ep^4*g*la*muu + 2*R^2*ep^4*g*muu^2)/(R^8*g^2*la^2 + 4*R^8*g^2*la*muu + 3*R^8*g^2*muu^2 + 2*R^8*g*la^2 + 8*R^8*g*la*muu + 10*R^8*g*muu^2 + R^8*la^2 + 4*R^8*la*muu + 3*R^8*muu^2 + 4*R^6*ep^2*g^2*la^2 + 8*R^6*ep^2*g^2*la*muu + 4*R^6*ep^2*g^2*muu^2 + 8*R^6*ep^2*g*la*muu + 8*R^6*ep^2*g*muu^2 - 4*R^6*ep^2*la^2 - 16*R^6*ep^2*la*muu - 12*R^6*ep^2*muu^2 - 6*R^4*ep^4*g^2*la^2 - 12*R^4*ep^4*g^2*la*muu - 6*R^4*ep^4*g^2*muu^2 - 12*R^4*ep^4*g*la*muu - 12*R^4*ep^4*g*muu^2 + 6*R^4*ep^4*la^2 + 24*R^4*ep^4*la*muu + 18*R^4*ep^4*muu^2 + 4*R^2*ep^6*g^2*la^2 + 12*R^2*ep^6*g^2*la*muu + 12*R^2*ep^6*g^2*muu^2 + 4*R^2*ep^6*g*la*muu - 4*R^2*ep^6*la^2 - 16*R^2*ep^6*la*muu - 12*R^2*ep^6*muu^2 + ep^8*g^2*la^2 + 4*ep^8*g^2*la*muu + 3*ep^8*g^2*muu^2 - 2*ep^8*g*la^2 - 8*ep^8*g*la*muu - 6*ep^8*g*muu^2 + ep^8*la^2 + 4*ep^8*la*muu + 3*ep^8*muu^2);
AR2 = (R^2*muu + R^2*g*la + R^2*g*muu)/(4*(R^2*muu - ep^2*muu + R^2*g*la + R^2*g*muu + ep^2*g*muu));
AR3 = -(R^4*ep^2*g^2*la^2 + 2*R^4*ep^2*g^2*la*muu + R^4*ep^2*g^2*muu^2 + 2*R^4*ep^2*g*la*muu + 2*R^4*ep^2*g*muu^2 - R^4*ep^2*la^2 - 4*R^4*ep^2*la*muu - 3*R^4*ep^2*muu^2 - R^2*ep^4*g^2*la^2 - 2*R^2*ep^4*g^2*la*muu - R^2*ep^4*g^2*muu^2 - 2*R^2*ep^4*g*la*muu - 2*R^2*ep^4*g*muu^2 + R^2*ep^4*la^2 + 4*R^2*ep^4*la*muu + 3*R^2*ep^4*muu^2)/(2*(R^8*g^2*la^2 + 4*R^8*g^2*la*muu + 3*R^8*g^2*muu^2 + 2*R^8*g*la^2 + 8*R^8*g*la*muu + 10*R^8*g*muu^2 + R^8*la^2 + 4*R^8*la*muu + 3*R^8*muu^2 + 4*R^6*ep^2*g^2*la^2 + 8*R^6*ep^2*g^2*la*muu + 4*R^6*ep^2*g^2*muu^2 + 8*R^6*ep^2*g*la*muu + 8*R^6*ep^2*g*muu^2 - 4*R^6*ep^2*la^2 - 16*R^6*ep^2*la*muu - 12*R^6*ep^2*muu^2 - 6*R^4*ep^4*g^2*la^2 - 12*R^4*ep^4*g^2*la*muu - 6*R^4*ep^4*g^2*muu^2 - 12*R^4*ep^4*g*la*muu - 12*R^4*ep^4*g*muu^2 + 6*R^4*ep^4*la^2 + 24*R^4*ep^4*la*muu + 18*R^4*ep^4*muu^2 + 4*R^2*ep^6*g^2*la^2 + 12*R^2*ep^6*g^2*la*muu + 12*R^2*ep^6*g^2*muu^2 + 4*R^2*ep^6*g*la*muu - 4*R^2*ep^6*la^2 - 16*R^2*ep^6*la*muu - 12*R^2*ep^6*muu^2 + ep^8*g^2*la^2 + 4*ep^8*g^2*la*muu + 3*ep^8*g^2*muu^2 - 2*ep^8*g*la^2 - 8*ep^8*g*la*muu - 6*ep^8*g*muu^2 + ep^8*la^2 + 4*ep^8*la*muu + 3*ep^8*muu^2));
AR4 = -(ep^2*(R^8*g^2*la^2 + 2*R^8*g^2*la*muu + R^8*g^2*muu^2 + 2*R^8*g*la*muu + 2*R^8*g*muu^2 - R^8*la^2 - 4*R^8*la*muu - 3*R^8*muu^2 + R^2*ep^6*g^2*la^2 + 4*R^2*ep^6*g^2*la*muu + 3*R^2*ep^6*g^2*muu^2 - 2*R^2*ep^6*g*la^2 - 8*R^2*ep^6*g*la*muu - 6*R^2*ep^6*g*muu^2 + R^2*ep^6*la^2 + 4*R^2*ep^6*la*muu + 3*R^2*ep^6*muu^2))/(2*(R^8*g^2*la^2 + 4*R^8*g^2*la*muu + 3*R^8*g^2*muu^2 + 2*R^8*g*la^2 + 8*R^8*g*la*muu + 10*R^8*g*muu^2 + R^8*la^2 + 4*R^8*la*muu + 3*R^8*muu^2 + 4*R^6*ep^2*g^2*la^2 + 8*R^6*ep^2*g^2*la*muu + 4*R^6*ep^2*g^2*muu^2 + 8*R^6*ep^2*g*la*muu + 8*R^6*ep^2*g*muu^2 - 4*R^6*ep^2*la^2 - 16*R^6*ep^2*la*muu - 12*R^6*ep^2*muu^2 - 6*R^4*ep^4*g^2*la^2 - 12*R^4*ep^4*g^2*la*muu - 6*R^4*ep^4*g^2*muu^2 - 12*R^4*ep^4*g*la*muu - 12*R^4*ep^4*g*muu^2 + 6*R^4*ep^4*la^2 + 24*R^4*ep^4*la*muu + 18*R^4*ep^4*muu^2 + 4*R^2*ep^6*g^2*la^2 + 12*R^2*ep^6*g^2*la*muu + 12*R^2*ep^6*g^2*muu^2 + 4*R^2*ep^6*g*la*muu - 4*R^2*ep^6*la^2 - 16*R^2*ep^6*la*muu - 12*R^2*ep^6*muu^2 + ep^8*g^2*la^2 + 4*ep^8*g^2*la*muu + 3*ep^8*g^2*muu^2 - 2*ep^8*g*la^2 - 8*ep^8*g*la*muu - 6*ep^8*g*muu^2 + ep^8*la^2 + 4*ep^8*la*muu + 3*ep^8*muu^2));
BR0 = -(R^2*(R^6*g^2*la^2 + 3*R^6*g^2*la*muu + 2*R^6*g^2*muu^2 + R^6*g*la^2 + 5*R^6*g*la*muu + 6*R^6*g*muu^2 - 3*R^2*ep^4*g^2*la^2 - 9*R^2*ep^4*g^2*la*muu - 6*R^2*ep^4*g^2*muu^2 + 3*R^2*ep^4*g*la^2 + 9*R^2*ep^4*g*la*muu + 6*R^2*ep^4*g*muu^2 + 4*ep^6*g^2*la^2 + 14*ep^6*g^2*la*muu + 12*ep^6*g^2*muu^2 - 4*ep^6*g*la^2 - 14*ep^6*g*la*muu - 12*ep^6*g*muu^2))/(R^8*g^2*la^2 + 4*R^8*g^2*la*muu + 3*R^8*g^2*muu^2 + 2*R^8*g*la^2 + 8*R^8*g*la*muu + 10*R^8*g*muu^2 + R^8*la^2 + 4*R^8*la*muu + 3*R^8*muu^2 + 4*R^6*ep^2*g^2*la^2 + 8*R^6*ep^2*g^2*la*muu + 4*R^6*ep^2*g^2*muu^2 + 8*R^6*ep^2*g*la*muu + 8*R^6*ep^2*g*muu^2 - 4*R^6*ep^2*la^2 - 16*R^6*ep^2*la*muu - 12*R^6*ep^2*muu^2 - 6*R^4*ep^4*g^2*la^2 - 12*R^4*ep^4*g^2*la*muu - 6*R^4*ep^4*g^2*muu^2 - 12*R^4*ep^4*g*la*muu - 12*R^4*ep^4*g*muu^2 + 6*R^4*ep^4*la^2 + 24*R^4*ep^4*la*muu + 18*R^4*ep^4*muu^2 + 4*R^2*ep^6*g^2*la^2 + 12*R^2*ep^6*g^2*la*muu + 12*R^2*ep^6*g^2*muu^2 + 4*R^2*ep^6*g*la*muu - 4*R^2*ep^6*la^2 - 16*R^2*ep^6*la*muu - 12*R^2*ep^6*muu^2 + ep^8*g^2*la^2 + 4*ep^8*g^2*la*muu + 3*ep^8*g^2*muu^2 - 2*ep^8*g*la^2 - 8*ep^8*g*la*muu - 6*ep^8*g*muu^2 + ep^8*la^2 + 4*ep^8*la*muu + 3*ep^8*muu^2);
BR2 = -(R^2*(R^6*g^2*la^2 + 4*R^6*g^2*la*muu + 3*R^6*g^2*muu^2 + 2*R^6*g*la^2 + 8*R^6*g*la*muu + 10*R^6*g*muu^2 + R^6*la^2 + 4*R^6*la*muu + 3*R^6*muu^2 - 3*R^2*ep^4*g^2*la^2 - 6*R^2*ep^4*g^2*la*muu - 3*R^2*ep^4*g^2*muu^2 - 6*R^2*ep^4*g*la*muu - 6*R^2*ep^4*g*muu^2 + 3*R^2*ep^4*la^2 + 12*R^2*ep^4*la*muu + 9*R^2*ep^4*muu^2 + 4*ep^6*g^2*la^2 + 12*ep^6*g^2*la*muu + 12*ep^6*g^2*muu^2 + 4*ep^6*g*la*muu - 4*ep^6*la^2 - 16*ep^6*la*muu - 12*ep^6*muu^2))/(2*(R^8*g^2*la^2 + 4*R^8*g^2*la*muu + 3*R^8*g^2*muu^2 + 2*R^8*g*la^2 + 8*R^8*g*la*muu + 10*R^8*g*muu^2 + R^8*la^2 + 4*R^8*la*muu + 3*R^8*muu^2 + 4*R^6*ep^2*g^2*la^2 + 8*R^6*ep^2*g^2*la*muu + 4*R^6*ep^2*g^2*muu^2 + 8*R^6*ep^2*g*la*muu + 8*R^6*ep^2*g*muu^2 - 4*R^6*ep^2*la^2 - 16*R^6*ep^2*la*muu - 12*R^6*ep^2*muu^2 - 6*R^4*ep^4*g^2*la^2 - 12*R^4*ep^4*g^2*la*muu - 6*R^4*ep^4*g^2*muu^2 - 12*R^4*ep^4*g*la*muu - 12*R^4*ep^4*g*muu^2 + 6*R^4*ep^4*la^2 + 24*R^4*ep^4*la*muu + 18*R^4*ep^4*muu^2 + 4*R^2*ep^6*g^2*la^2 + 12*R^2*ep^6*g^2*la*muu + 12*R^2*ep^6*g^2*muu^2 + 4*R^2*ep^6*g*la*muu - 4*R^2*ep^6*la^2 - 16*R^2*ep^6*la*muu - 12*R^2*ep^6*muu^2 + ep^8*g^2*la^2 + 4*ep^8*g^2*la*muu + 3*ep^8*g^2*muu^2 - 2*ep^8*g*la^2 - 8*ep^8*g*la*muu - 6*ep^8*g*muu^2 + ep^8*la^2 + 4*ep^8*la*muu + 3*ep^8*muu^2));
BR3 = -(R^2*(ep^2*muu - ep^2*g*muu))/(2*(R^2*muu - ep^2*muu + R^2*g*la + R^2*g*muu + ep^2*g*muu));
BR4 = -(R^4*(g - 1)*(ep^8*g*la^2 - 3*ep^8*muu^2 - ep^8*la^2 + 3*ep^8*g*muu^2 + R^4*ep^4*la^2 + 3*R^4*ep^4*muu^2 - 4*ep^8*la*muu + R^4*ep^4*g*la^2 + R^4*ep^4*g*muu^2 + 4*ep^8*g*la*muu + 4*R^4*ep^4*la*muu + 2*R^4*ep^4*g*la*muu))/(2*(R^8*g^2*la^2 + 4*R^8*g^2*la*muu + 3*R^8*g^2*muu^2 + 2*R^8*g*la^2 + 8*R^8*g*la*muu + 10*R^8*g*muu^2 + R^8*la^2 + 4*R^8*la*muu + 3*R^8*muu^2 + 4*R^6*ep^2*g^2*la^2 + 8*R^6*ep^2*g^2*la*muu + 4*R^6*ep^2*g^2*muu^2 + 8*R^6*ep^2*g*la*muu + 8*R^6*ep^2*g*muu^2 - 4*R^6*ep^2*la^2 - 16*R^6*ep^2*la*muu - 12*R^6*ep^2*muu^2 - 6*R^4*ep^4*g^2*la^2 - 12*R^4*ep^4*g^2*la*muu - 6*R^4*ep^4*g^2*muu^2 - 12*R^4*ep^4*g*la*muu - 12*R^4*ep^4*g*muu^2 + 6*R^4*ep^4*la^2 + 24*R^4*ep^4*la*muu + 18*R^4*ep^4*muu^2 + 4*R^2*ep^6*g^2*la^2 + 12*R^2*ep^6*g^2*la*muu + 12*R^2*ep^6*g^2*muu^2 + 4*R^2*ep^6*g*la*muu - 4*R^2*ep^6*la^2 - 16*R^2*ep^6*la*muu - 12*R^2*ep^6*muu^2 + ep^8*g^2*la^2 + 4*ep^8*g^2*la*muu + 3*ep^8*g^2*muu^2 - 2*ep^8*g*la^2 - 8*ep^8*g*la*muu - 6*ep^8*g*muu^2 + ep^8*la^2 + 4*ep^8*la*muu + 3*ep^8*muu^2));
 


r = 0:0.02:ep;
t = 0:pi/120:2*pi;
[r,t] = meshgrid(r,t);

ur1 = -al*(S*r.*(2*AR1*la*cos(2*t).*r.^2 - 2*AR0*muu + BR0*la*cos(2*t) + BR0*muu*cos(2*t)))/(2*g*muu*(la + muu));
ut1 = al*(S*r.*sin(2*t).*(BR0*la + BR0*muu + 4*AR1*la*r.^2 + 6*AR1*muu*r.^2))/(2*g*muu*(la + muu));

ux1 = cos(t).*ur1 - sin(t).*ut1;
uy1 = sin(t).*ur1 + cos(t).*ut1;

Srr1 = S*(2*AR0 - BR0*cos(2*t));
Srt1 = S*sin(2*t).*(6*AR1*r.^2 + BR0);
Stt1 = S*(12*AR1*cos(2*t).*r.^2 + 2*AR0 + BR0*cos(2*t));

Sxx1 = (Srr1.*cos(t) - Srt1.*sin(t)).*cos(t) + (Srt1.*cos(t) - Stt1.*sin(t)).*(-sin(t));
Sxy1 = (Srr1.*cos(t) - Srt1.*sin(t)).*sin(t) + (Srt1.*cos(t) - Stt1.*sin(t)).*cos(t);
Syy1 = (Srr1.*sin(t) + Srt1.*cos(t)).*sin(t) + (Srt1.*sin(t) + Stt1.*cos(t)).*cos(t);

x = r.*cos(t);
y = r.*sin(t);
ux = figure;
surf(x,y,ux1);
title('ux');
uy = figure;
surf(x,y,uy1);
title('uy');
Sxx = figure;
surf(x,y,Sxx1);
title('Sxx');
Sxy = figure;
surf(x,y,Sxy1);
title('Sxy');
Syy = figure;
surf(x,y,Syy1);
title('Syy');



r = ep:0.02:R;
t = 0:pi/120:2*pi;
[r,t] = meshgrid(r,t);

ur2 = -al*(S*(BR4*la*cos(2*t) + BR4*muu*cos(2*t) + BR3*la*r.^2 - 2*AR2*muu*r.^4 + BR3*muu*r.^2 - 2*AR4*la*r.^2.*cos(2*t) + 2*AR3*la*r.^6.*cos(2*t) + BR2*la*r.^4.*cos(2*t) - 4*AR4*muu*r.^2.*cos(2*t) + BR2*muu*r.^4.*cos(2*t)))./(2*muu*r.^3*(la + muu));
ut2 = al*(S*sin(2*t).*(4*AR3*la*r.^6 - BR4*muu - BR4*la + BR2*la*r.^4 - 2*AR4*muu*r.^2 + 6*AR3*muu*r.^6 + BR2*muu*r.^4))./(2*muu*r.^3*(la + muu));

ux2 = cos(t).*ur2 - sin(t).*ut2;
uy2 = sin(t).*ur2 + cos(t).*ut2;

Srr2 = (S*(2*AR2*r.^4 + BR3*r.^2 + 3*BR4*cos(2*t) - 4*AR4*r.^2.*cos(2*t) - BR2*r.^4.*cos(2*t)))./r.^4;
Srt2 = (S*sin(2*t).*(6*AR3*r.^6 + BR2*r.^4 - 2*AR4*r.^2 + 3*BR4))./r.^4;
Stt2 = (S*(2*AR2*r.^4 - BR3*r.^2 - 3*BR4*cos(2*t) + 12*AR3*r.^6.*cos(2*t) + BR2*r.^4.*cos(2*t)))./r.^4;

Sxx2 = (Srr2.*cos(t) - Srt2.*sin(t)).*cos(t) + (Srt2.*cos(t) - Stt2.*sin(t)).*(-sin(t));
Sxy2 = (Srr2.*cos(t) - Srt2.*sin(t)).*sin(t) + (Srt2.*cos(t) - Stt2.*sin(t)).*cos(t);
Syy2 = (Srr2.*sin(t) + Srt2.*cos(t)).*sin(t) + (Srt2.*sin(t) + Stt2.*cos(t)).*cos(t);

x = r.*cos(t);
y = r.*sin(t);
figure(ux);
hold on
surf(x,y,ux2);
figure(uy);
hold on
surf(x,y,uy2);
figure(Sxx);
hold on
surf(x,y,Sxx2);
figure(Sxy);
hold on
surf(x,y,Sxy2);
figure(Syy);
hold on
surf(x,y,Syy2);


%% Resultados XFEM


X=ES.Mnodo(:,2);%Coord X de los nodos

Y=ES.Mnodo(:,3); %Coord Y de los nodos

U=ES.U(1:2:(2*ES.Nnodo-1)); % Despla X de los nodos

V=ES.U(2:2:(2*ES.Nnodo)); % Despla Y de los nodos

Psi = ES.psi; % Level set de los nodos

Ux=0*U; % Exten x de los nodos (Los -1e13 son donde no hay)

Vx=0*V; % Exten y de los nodos
for n=1:ES.Nnodo
    ne=ES.NGLX(n);
    if ne==0
        Ux(n)=1e13;
        Vx(n)=1e13;
    else
        Ux(n)=ES.U(ne);
        Vx(n)=ES.U(ne+1);
    end
end

%% Plot de desplazamietnos



paso=2e-3;
vecX=-.2:paso:1.999e-1;
vecY=-2:paso:1.999;

xplot=zeros(length(vecX)*length(vecY),1);
yplot=xplot;
ufem=xplot;
vfem=yplot;
uxfem=xplot;
vxfem=yplot;
UTot=xplot;
VTot=yplot;
uanal = xplot; vanal = yplot;

SigmaX_xfem=xplot;
SigmaY_xfem=xplot;
SigmaXY_xfem=xplot;
SigmaXanal = xplot;
SigmaYanal = xplot;
SigmaXYanal = xplot;


xm = min(ES.Mnodo(:,2));
ym = min(ES.Mnodo(:,3));
dx = ES.Lx/ES.Nx;
dy = ES.Ly/ES.Ny;

punto=1;

for i = 1:length(vecX)
    for j = 1:length(vecY)
        xplot(punto)=vecX(i);
        yplot(punto)=vecY(j);
        
        indx = 1+floor( (vecX(i)-xm)/dx);
        indy = 1+floor( (vecY(j)-ym)/dy);
        
        xizq = xm + dx*(indx-1);
        yaba = ym + dy*(indy-1);
        
        if vecY(j)-yaba <= vecX(i) - xizq %Elemento abajo a la izquierda
            ind1=indx + (indy-1)*(ES.Nx+1);
            ind2= ind1+1;
            ind3 = ind2 + ES.Nx+1;
            xi = (vecY(j)-yaba)/dy;
            eta = (vecX(i)-xizq)/dx-xi;
        else
            ind1= indx + (indy-1)*(ES.Nx+1);
            ind3=ind1+ES.Nx+1;
            ind2=ind3+1;
            eta=(vecX(i)-xizq)/dx;
            xi=(vecY(j)-yaba)/dy - eta;
        end
        
       
        
        uele = [U(ind1),U(ind2),U(ind3)]';
        vele = [V(ind1),V(ind2),V(ind3)]';
        uxele = [Ux(ind1),Ux(ind2),Ux(ind3)]';
        vxele = [Vx(ind1),Vx(ind2),Vx(ind3)]';
        
        N1 = 1-eta-xi;
        N2 = eta;
        N3 = xi;
        N=[N1 N2 N3];
        
        psiele = [Psi(ind1),Psi(ind2),Psi(ind3)]';
        apsi = abs(psiele);
        gxfem = N*apsi - abs(N*psiele);
        w = N*gxfem;
        
        ufem(punto) = N*uele;% + w *uxele;
        vfem(punto) = N*vele;%;
        
        uxfem(punto)= w*uxele;
        vxfem(punto)= w*vxele;
        
        UTot(punto)=ufem(punto)+uxfem(punto);
        VTot(punto)=vfem(punto)+vxfem(punto);
        
        x=[X(ind1);X(ind2);X(ind3)];
        y = [Y(ind1);Y(ind2);Y(ind3)];
        dN_detachi = [-1 1 0; -1 0 1];
        J=dN_detachi*[x,y];
        
        dN_dxy=J\dN_detachi;
        
        if N*psiele >0
            signo=1;
        else
            signo=-1;
        end
        dg_detachi = dN_detachi * apsi - signo*dN_detachi*psiele;
        
        dw_detachi = gxfem * dN_detachi + dg_detachi*N;
        dw_dxy = J\dw_detachi;
        
        B(1,1:2:5) = dN_dxy(1,:);
        B(2,2:2:6) = dN_dxy(2,:);
        B(3,1:2:5) = dN_dxy(2,:);
        B(3,2:2:6) = dN_dxy(1,:);
        B(1,7:2:11) = dw_dxy(1,:);
        B(2,8:2:12) = dw_dxy(2,:);
        B(3,7:2:11) = dw_dxy(2,:);
        B(3,8:2:12) = dw_dxy(1,:);
        
        uXfemizado = [uele(1);vele(1);uele(2);vele(2);uele(3);vele(3);
                      uxele(1);vxele(1);uxele(2);vxele(2);uxele(3);vxele(3)];
        epsss = B*uXfemizado;
        epsX = epsss(1); epsY=epsss(2); epsXY = epsss(3);
        r = sqrt(vecX(i).^2 + vecY(j).^2);
         t = atan2(vecY(j),vecX(i)); if t<0; t=t+2*pi; end
         
         if r<=ep
             ur1 = -al*(S*r.*(2*AR1*la*cos(2*t).*r.^2 - 2*AR0*muu + BR0*la*cos(2*t) + BR0*muu*cos(2*t)))/(2*g*muu*(la + muu));
ut1 = al*(S*r.*sin(2*t).*(BR0*la + BR0*muu + 4*AR1*la*r.^2 + 6*AR1*muu*r.^2))/(2*g*muu*(la + muu));

uanal(punto) = cos(t).*ur1 - sin(t).*ut1;
vanal(punto) = sin(t).*ur1 + cos(t).*ut1;
Srr1 = S*(2*AR0 - BR0*cos(2*t));
Srt1 = S*sin(2*t).*(6*AR1*r.^2 + BR0);
Stt1 = S*(12*AR1*cos(2*t).*r.^2 + 2*AR0 + BR0*cos(2*t));


             SigmaX_xfem(punto) = (la + 2*muu) * epsX + la*epsY;
             SigmaY_xfem(punto) = la*epsX + (la+2*muu)*epsY;
             SigmaXY_xfem(punto) = muu*epsXY;

SigmaXanal(punto) = (Srr1.*cos(t) - Srt1.*sin(t)).*cos(t) + (Srt1.*cos(t) - Stt1.*sin(t)).*(-sin(t));
SigmaXYanal(punto) = (Srr1.*cos(t) - Srt1.*sin(t)).*sin(t) + (Srt1.*cos(t) - Stt1.*sin(t)).*cos(t);
SigmaYanal(punto) = (Srr1.*sin(t) + Srt1.*cos(t)).*sin(t) + (Srt1.*sin(t) + Stt1.*cos(t)).*cos(t);

         else
             
    ur2 = -al*(S*(BR4*la*cos(2*t) + BR4*muu*cos(2*t) + BR3*la*r.^2 - 2*AR2*muu*r.^4 + BR3*muu*r.^2 - 2*AR4*la*r.^2.*cos(2*t) + 2*AR3*la*r.^6.*cos(2*t) + BR2*la*r.^4.*cos(2*t) - 4*AR4*muu*r.^2.*cos(2*t) + BR2*muu*r.^4.*cos(2*t)))./(2*muu*r.^3*(la + muu));
ut2 = al*(S*sin(2*t).*(4*AR3*la*r.^6 - BR4*muu - BR4*la + BR2*la*r.^4 - 2*AR4*muu*r.^2 + 6*AR3*muu*r.^6 + BR2*muu*r.^4))./(2*muu*r.^3*(la + muu));

uanal(punto) = cos(t).*ur2 - sin(t).*ut2;
vanal(punto) = sin(t).*ur2 + cos(t).*ut2;


             SigmaX_xfem(punto) = ( (la + 2*muu) * epsX + la*epsY ) /al;
             SigmaY_xfem(punto) = (la*epsX + (la+2*muu)*epsY ) /al;
             SigmaXY_xfem(punto) = (muu*epsXY ) /al;

             
             Srr2 = (S*(2*AR2*r.^4 + BR3*r.^2 + 3*BR4*cos(2*t) - 4*AR4*r.^2.*cos(2*t) - BR2*r.^4.*cos(2*t)))./r.^4;
Srt2 = (S*sin(2*t).*(6*AR3*r.^6 + BR2*r.^4 - 2*AR4*r.^2 + 3*BR4))./r.^4;
Stt2 = (S*(2*AR2*r.^4 - BR3*r.^2 - 3*BR4*cos(2*t) + 12*AR3*r.^6.*cos(2*t) + BR2*r.^4.*cos(2*t)))./r.^4;

SigmaXanal(punto) = (Srr2.*cos(t) - Srt2.*sin(t)).*cos(t) + (Srt2.*cos(t) - Stt2.*sin(t)).*(-sin(t));
SigmaXYanal(punto)= (Srr2.*cos(t) - Srt2.*sin(t)).*sin(t) + (Srt2.*cos(t) - Stt2.*sin(t)).*cos(t);
SigmaYanal(punto) = (Srr2.*sin(t) + Srt2.*cos(t)).*sin(t) + (Srt2.*sin(t) + Stt2.*cos(t)).*cos(t);

         end
        
        
        punto=punto+1;
    end
end


figure
scatter3(xplot,yplot,UTot,100,UTot,'filled')
title('UXFEM')
figure
scatter3(xplot,yplot,VTot,10,VTot)
title('VXFEM')
figure
scatter3(xplot,yplot,uanal,10,uanal,'filled')
title('U analitico')
figure
scatter3(xplot,yplot,vanal,10,vanal,'filled');
title('V analitico')
ErrorURelat = 0*uanal;
ErrorVRelat= 0*vanal;
IndiceU = 0.0001*mean(abs(uanal));
IndiceV=0.0001*mean(abs(vanal));
for punto = 1:length(ErrorURelat)
    if abs(uanal(punto))>IndiceU
        ErrorURelat(punto)=1-UTot(punto)/uanal(punto);
    else
        ErrorURelat(punto)=NaN;
    end
    if abs(vanal(punto))>IndiceV
        ErrorVRelat(punto)=1-VTot(punto)/vanal(punto);
    else
        ErrorVRelat(punto)=NaN;
    end
end
figure
scatter3(xplot,yplot,ErrorURelat,10,ErrorURelat,'filled');
title('Error U - Relat')
figure
scatter3(xplot,yplot,ErrorVRelat,10,ErrorVRelat,'filled');
title('Error V - Relat')

figure
scatter3(xplot,yplot,SigmaX_xfem,100,SigmaX_xfem,'filled')
title('\sigma_x xfem')

figure
scatter3(xplot,yplot,SigmaY_xfem,100,SigmaY_xfem,'filled')
title('\sigma_y xfem')


figure
scatter3(xplot,yplot,SigmaXY_xfem,100,SigmaXY_xfem,'filled')
title('\tau_{xy} xfem')





figure
scatter3(xplot,yplot,SigmaXanal,100,SigmaXanal,'filled')
title('\sigma_x anal')

figure
scatter3(xplot,yplot,SigmaYanal,100,SigmaYanal,'filled')
title('\sigma_y anal')


figure
scatter3(xplot,yplot,SigmaXYanal,100,SigmaXYanal,'filled')
title('\tau_{xy} anal')



ErrorSXRelat = 0*SigmaXYanal;
ErrorSYRelat= 0*SigmaXYanal;
ErrorTXYRelat = 0*SigmaXYanal;
IndiceSX = 0.0001*mean(abs(SigmaXanal));
IndiceSY=0.0001*mean(abs(SigmaYanal));
IndiceTXY = 0.0001*mean(abs(SigmaXYanal));
for punto = 1:length(ErrorSXRelat)
    if abs(SigmaXanal(punto))>IndiceSX
        ErrorSXRelat(punto)=1-SigmaX_xfem(punto)/SigmaXanal(punto);
    else
        ErrorSXRelat(punto)=NaN;
    end
    if abs(SigmaYanal(punto))>IndiceSY
        ErrorSYRelat(punto)=1-SigmaY_xfem(punto)/SigmaYanal(punto);
    else
        ErrorSYRelat(punto)=NaN;
    end
    if abs(SigmaXYanal(punto))>IndiceTXY
        ErrorTXYRelat(punto)=1-SigmaXY_xfem(punto)/SigmaXYanal(punto);
    else
        ErrorTXYRelat(punto)=NaN;
    end
end
figure
scatter3(xplot,yplot,ErrorSXRelat,10,ErrorSXRelat,'filled');
title('Error \sigma_x - Relat')


figure
scatter3(xplot,yplot,ErrorSYRelat,10,ErrorSYRelat,'filled');
title('Error \sigma_y - Relat')

figure
scatter3(xplot,yplot,ErrorTXYRelat,10,ErrorTXYRelat,'filled');
title('Error \tau - Relat')


return

%% Tensiones

xplot=zeros((length(vecX)-1)*(length(vecY)-1),1);
yplot=xplot;
punto=1;

xvec=vecX; yvec=vecY;



for i=1:length(xvec)-1
    for j=1:length(yvec)-1
        
        xplot(punto)=xvec(i+1)*0.5+xvec(i)*0.5;
        yplot(punto) = yvec(j+1)*0.5+yvec(j)*0.5;
        dx = xvec(i+1)-xvec(i);
        dy = yvec(j+1)-yvec(j);
        
        deltaX = 0.5*UTot( j+1 +(i+1-1)*length(yvec)) + 0.5*UTot( j +(i+1-1)*length(yvec) ) - 0.5*UTot( j+1 + (i-1)*length(yvec)) - 0.5*UTot( j + (i-1)*length(yvec));
        epsX = deltaX/dx;
        
        deltaY = 0.5*VTot( j+1 + (i+1-1)*length(yvec)) + 0.5*VTot( j+1 + (i-1)*length(yvec) ) - 0.5*VTot( j + (i+1-1)*length(yvec) ) - 0.5*VTot(j + (i-1)*length(yvec));
        epsY = deltaY/dy;
        
        deltaX = 0.5*UTot( j+1 + (i+1-1)*length(yvec)) + 0.5*UTot( j+1 + (i-1)*length(yvec) ) - 0.5*UTot( j + (i+1-1)*length(yvec) ) - 0.5*UTot(j + (i-1)*length(yvec));
        deltaY = 0.5*VTot( j+1 +(i+1-1)*length(yvec)) + 0.5*VTot( j +(i+1-1)*length(yvec) ) - 0.5*VTot( j+1 + (i-1)*length(yvec)) - 0.5*VTot( j + (i-1)*length(yvec));
        epsXY = deltaY/dx + deltaX/dy;
        
        
        r = sqrt(xplot(punto).^2 + yplot(punto).^2);
         t = atan2(yplot(punto),xplot(punto)); if t<0; t=t+2*pi; end
        
        if r <=ep 
             
             
             
        else
             
        end
        
        
        punto=punto+1;
    end
end

