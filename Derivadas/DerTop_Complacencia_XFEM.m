function ES = DerTop_Complacencia_XFEM(ES)

% Rutina para hallar la derivada topologica de la complacencia
% segun el paper de Lopes 2015. Recibe un x2 respecto a la del paper pues
% la del paper es La energia potencial total que es -0.5 Complacencia.
%
% Calcuando las derivadas especiales en los elementos finitos y luego
% asignando al nodo un valor ponderado (VP) segun la siguiente formula:
%
% VP = 6 * \int N * DT.
% Siendo N la funcion de forma del nodo en el elemento, DT la derivada
% topologica y la integral en el elemento de referencia ( (0,0) (1,0) (0,1) )
% El valor 6 es porque si DT fuera cte = 1 quisieramos que VP = 1 y la
% integral de N *1 = 1/6
%
% Como algunos (la mayoria) de los nodos tienen mas de un elemento, se
% calcula el VP de cada elemento. Luego para definir la DT del
% nodo se calcula la media de los VP.
%
% Devuelve o actualiza:
%
% ES.DTC:     Derivada topologica de la complacencia. Tiene tanta entradas como nodo la malla

% =========================================================================
% === Primeras variables ==================================================
% =========================================================================

SumVP = zeros(ES.Nnodo,1); 
EleCounter = zeros(ES.Nnodo,1);

% =========================================================================
% === Iteracion en los elementos ==========================================
% =========================================================================

for ele = 1:ES.Nelem
    
    ne = ES.Melem(ele,3:5);

    Xe = ES.Mnodo(ne,2);
    Ye = ES.Mnodo(ne,3);

    dN_detachi=[-1 1 0; -1 0 1];

    J = dN_detachi*[Xe,Ye];

    dN_dxy = J\dN_detachi; 

    if ES.EI(ele)
       
        psie = ES.psi(ne); 
        
        ArX = ES.EGLX(ele,:); 
        NArX = ES.AriX(ArX,:); 
        
        Uele = zeros(6,1); 
        Uele(1:2:5) = ES.U(2*ne-1);
        Uele(2:2:6) = ES.U(2*ne); 
		
        Uxele = zeros(12,1); 
        ngx = full(ES.NodAriX(ne,ArX));
        Uxele(1:2:11) = ES.U(ngx);
        Uxele(2:2:12) = ES.U(ngx+1); 

        UGLele = [Uele;Uxele];
		
        if psie(1)*psie(2)<0 
            etaR1 = psie(1)/(psie(1)-psie(2));

            if ( ne(1)==NArX(1,1) || ne(1)==NArX(1,2) ) && ( ne(2)==NArX(1,1) || ne(2) == NArX(1,2) ) 
                Ar1R1 = 1; 
                Ar2R1 = 0; 
            else
                Ar2R1 = 1; 
                Ar1R1 = 0;
            end
        else
            etaR1 = 0.5;
            Ar1R1 = 0; Ar2R1 = 0;
        end
		
        if psie(2)*psie(3)<0
            etaR2 = psie(3)/(psie(3)-psie(2));
            if ( ne(2)==NArX(1,1) || ne(2)==NArX(1,2) ) && ( ne(3)==NArX(1,1) || ne(3) == NArX(1,2) ) 
                Ar1R2 = 1;
                Ar2R2 = 0; 
            else
                Ar2R2 = 1; 
                Ar1R2 = 0;
            end
        else
            etaR2 = 0.5;
            Ar1R2 = 0; Ar2R2 = 0;
        end  
    
        if psie(3)*psie(1)<0
            xiR3 = psie(1)/(psie(1)-psie(3));
            if ( ne(1)==NArX(1,1) || ne(1)==NArX(1,2) ) && ( ne(3)==NArX(1,1) || ne(3) == NArX(1,2) )
                Ar1R3 = 1;
                Ar2R3 = 0; 
            else
                Ar2R3 = 1; 
                Ar1R3 = 0;
            end
        else
            xiR3 = 0.5;
            Ar1R3 = 0; Ar2R3 = 0;
        end  

        etaG = [1/3 0.6 0.2 0.2];
        xiG = [1/3 0.2 0.6 0.2];
        PesoG = [-0.28125 0.260416666666666667 0.260416666666666667  0.260416666666666667];

        etaSubElements = [0      1      0   etaR1;
                         etaR1 etaR2   0   etaR2;
                          0     etaR1 etaR2   0  ];
   
        xiSubElements = [0      0       1      0    ;
                         0   1-etaR2   xiR3  1-etaR2;
                        xiR3   0    1-etaR2  xiR3  ];

        
        Ar1SE = [0      0      0     Ar1R1;
                 Ar1R1  Ar1R2  Ar1R3 Ar1R2;
                 Ar1R3  Ar1R1  Ar1R2 Ar1R3 ]; 
    
        Ar2SE = [0      0      0     Ar2R1;
                 Ar2R1  Ar2R2  Ar2R3 Ar2R2;
                 Ar2R3  Ar2R1  Ar2R2 Ar2R3 ]; 
				 
        for subEl = 1:4 
		
            etase = etaSubElements(:,subEl); 
            xise = xiSubElements(:,subEl);
            Ar1 = Ar1SE(:,subEl); 
            Ar2 = Ar2SE(:,subEl); 
            
            JacSe = dN_detachi*[etase xise]; 
            
            Jse = det( JacSe ); 
			
            dNse_detachi = JacSe\dN_detachi; 
            
			
            etaeval =  [1/3 1/3 1/3]* etase;
            xieval = [1/3 1/3 1/3] *xise; 
            
            psiRep = [1-etaeval-xieval,etaeval,xieval]*psie; 
            
            if psiRep<0
                E=ES.PropMat(1,2); 
                nu=ES.PropMat(1,3); 
                gammaP = 1/ES.gamma ;
                b=ES.CB.Neu.Vol(1,:);  
            else
                E=ES.PropMat(2,2);
                nu=ES.PropMat(2,3); 
                gammaP = ES.gamma ;
                b=ES.CB.Neu.Vol(2,:); 
            end

            if ES.PEL==1
               mu = E/(2*(1+nu));
               lambda = nu*E/(1-nu^2);
            elseif ES.PEL==2
               mu = E/(2*(1+nu));
               lambda = nu*E/( ( 1+nu)*(1-2*nu));
               E = E/(1-nu^2);
               nu = nu / (1-nu);
            end

            C = (E/(1-nu^2))* [ 1 nu 0;
                            nu 1 0;
                            0 0 0.5*(1-nu)];

            alpha1= (lambda+mu)/mu;
            alpha2= (lambda+3*mu)/(lambda+mu);
            
            EleCounter(ne) = EleCounter(ne) + Jse.*ones(3,1);

            for puntoG=1:4 
			
                
                N1se = 1-etaG(puntoG) - xiG(puntoG); 
                N2se = etaG(puntoG); 
                N3se = xiG(puntoG); 
                Nse = [N1se,N2se,N3se];
                etaeval = Nse * etase; 
                xieval = Nse *xise; 
				
                B = zeros(3,18);
                
				
                B(1,1:2:5)=dN_dxy(1,:);
                B(2,2:2:6)=dN_dxy(2,:);
                B(3,2:2:6)=dN_dxy(1,:);
                B(3,1:2:5)=dN_dxy(2,:);
                                           
										   
                N1 = 1-etaeval-xieval; 
                N2 = etaeval; 
                N3 = xieval; 
                N = [N1,N2,N3]; 
				
                g1 = Nse * Ar1;
                g2 = Nse * Ar2;
                
				
                w1 = N*g1;
                w2 = N*g2;
				
                dg1_detachi = dNse_detachi*Ar1; 
                dg2_detachi = dNse_detachi*Ar2;
                dw1_detachi = dg1_detachi*N + dN_detachi*g1;
                dw2_detachi = dg2_detachi*N + dN_detachi*g2;
                dw1_dxy = J\dw1_detachi; 
                dw2_dxy = J\dw2_detachi;
                
				
                B(1,7:2:11)=dw1_dxy(1,:);
                B(2,8:2:12)=dw1_dxy(2,:);
                B(3,8:2:12)=dw1_dxy(1,:);
                B(3,7:2:11)=dw1_dxy(2,:);
                
                B(1,13:2:17) = dw2_dxy(1,:);
                B(2,14:2:18) = dw2_dxy(2,:);
                B(3,14:2:18) = dw2_dxy(1,:);
                B(3,13:2:17) = dw2_dxy(2,:);


                Eps = B *UGLele;

                Sigma=C*Eps;
				
                Sigma=[Sigma(1), Sigma(3); Sigma(3), Sigma(2)];

                trS = Sigma(1,1) + Sigma(2,2);

                SigmaP = 0.5*( (1-gammaP)/(1+gammaP*alpha2) ) * ( (1+alpha2)*Sigma + 0.5*(alpha1-alpha2)*( (1-gammaP)/(1+gammaP*alpha1) )*trS*eye(2) );

                SigmaP = [SigmaP(1,1);SigmaP(2,2);SigmaP(1,2)];

                Ux = [ N , w1 , w2 ] * UGLele(1:2:17) ; 
                Uy =  [ N , w1 , w2 ] * UGLele(2:2:18) ;

                DerTopG = SigmaP'*Eps + (1-gammaP) * (b(1) * Ux + b(2)*Uy );
              
                DerTopG = 2 * DerTopG;

                ValPerNodo =  ( N * DerTopG )'; 

                SumVP(ne) = SumVP(ne)+ 6 * Jse *PesoG(puntoG) * ValPerNodo;



            end


        end

    elseif ES.VA(ele)

        EleCounter(ne) = EleCounter(ne) + ones(3,1);
		
        E = ES.PropMat(1,2); 
        nu = ES.PropMat(1,3); 

        if ES.PEL==1
           mu = E/(2*(1+nu));
           lambda = nu*E/(1-nu^2);
        elseif ES.PEL==2
           mu = E/(2*(1+nu));
           lambda = nu*E/( ( 1+nu)*(1-2*nu));
           E = E/(1-nu^2);
           nu = nu / (1-nu);
        end

        Uele = zeros(6,1); 
        Uele(1:2:5) = ES.U(2*ne-1); 
        Uele(2:2:6) = ES.U(2*ne); 

        B = zeros(3,6);
        B(1,1:2:5)=dN_dxy(1,:);
        B(2,2:2:6)=dN_dxy(2,:);
        B(3,2:2:6)=dN_dxy(1,:);
        B(3,1:2:5)=dN_dxy(2,:);

        Eps = B*Uele;

        C = (E/(1-nu^2))* [ 1 nu 0;
                        nu 1 0;
                        0 0 0.5*(1-nu)];

        Sigma=C*Eps;
		
        Sigma=[Sigma(1), Sigma(3); Sigma(3), Sigma(2)];

        trS = Sigma(1,1) + Sigma(2,2); 

        gammaP = 1/ES.gamma ; 
		
        alpha1= (lambda+mu)/mu;
        alpha2= (lambda+3*mu)/(lambda+mu);

        SigmaP = 0.5*( (1-gammaP)/(1+gammaP*alpha2) ) * ( (1+alpha2)*Sigma + 0.5*(alpha1-alpha2)*( (1-gammaP)/(1+gammaP*alpha1) )*trS*eye(2) );

        SigmaP = [SigmaP(1,1);SigmaP(2,2);SigmaP(1,2)];

        b=ES.CB.Neu.Vol(1,:);


        etaG = [0.5 0 0.5]; 
        xiG = [0.5 0.5 0];
        PesoG = [1/6 1/6 1/6];


        for puntoG = 1:3 
		
            N1 = 1-etaG(puntoG) -xiG(puntoG) ;
            N2 = etaG(puntoG); 
            N3 = xiG(puntoG); 
            N=[N1,N2,N3]; 

            Ux = N * Uele(1:2:5) ; 
            Uy = N * Uele(2:2:6) ; 

            DerTopG = dot(SigmaP, Eps) + (1-gammaP) * (b(1) * Ux + b(2)*Uy );
          
            DerTopG = 2 * DerTopG;

            ValPerNodo = ( N * DerTopG )'; 

            SumVP(ne) =SumVP(ne)+ 6 * PesoG(puntoG) * ValPerNodo; 

        end 

    else

        EleCounter(ne) = EleCounter(ne) + ones(3,1);

        E = ES.PropMat(2,2); 
        nu = ES.PropMat(2,3); 

        if ES.PEL==1
           mu = E/(2*(1+nu));
           lambda = nu*E/(1-nu^2);
        elseif ES.PEL==2
           mu = E/(2*(1+nu));
           lambda = nu*E/( ( 1+nu)*(1-2*nu));
           E = E/(1-nu^2);
           nu = nu / (1-nu);
        end

        Uele = zeros(6,1);
        Uele(1:2:5) = ES.U(2*ne-1); 
        Uele(2:2:6) = ES.U(2*ne); 

        B = zeros(3,6);
        B(1,1:2:5)=dN_dxy(1,:);
        B(2,2:2:6)=dN_dxy(2,:);
        B(3,2:2:6)=dN_dxy(1,:);
        B(3,1:2:5)=dN_dxy(2,:);

        Eps = B*Uele;


        C = (E/(1-nu^2))* [ 1 nu 0;
                        nu 1 0;
                        0 0 0.5*(1-nu)];
        Sigma=C*Eps;

        Sigma=[Sigma(1), Sigma(3); Sigma(3), Sigma(2)];

        trS = Sigma(1,1) + Sigma(2,2); 

        gammaP = ES.gamma ; 
		
        alpha1= (lambda+mu)/mu;
        alpha2= (lambda+3*mu)/(lambda+mu);
		
        SigmaP = 0.5*( (1-gammaP)/(1+gammaP*alpha2) ) * ( (1+alpha2)*Sigma + 0.5*(alpha1-alpha2)*( (1-gammaP)/(1+gammaP*alpha1) )*trS*eye(2) );


        SigmaP = [SigmaP(1,1);SigmaP(2,2);SigmaP(1,2)];

        b=ES.CB.Neu.Vol(2,:);

        etaG = [0.5 0 0.5];
        xiG = [0.5 0.5 0]; 
        PesoG = [1/6 1/6 1/6];
		

        for puntoG = 1:3 

            N1 = 1-etaG(puntoG) -xiG(puntoG) ;
            N2 = etaG(puntoG); 
            N3 = xiG(puntoG); 
            N=[N1,N2,N3]; 
			
            Ux = N  * Uele(1:2:5) ; 
            Uy = N  * Uele(2:2:6) ; 

            DerTopG = dot(SigmaP, Eps) + (1-gammaP) * (b(1) * Ux + b(2)*Uy ); 
            DerTopG = 2 * DerTopG;
            
            ValPerNodo = ( N * DerTopG )'; 

            SumVP(ne) = SumVP(ne)+6 * PesoG(puntoG) * ValPerNodo;

        end 

    end

end

ES.DTC = SumVP./EleCounter;
