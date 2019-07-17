
disp('SLM-SLP with adaptive scheme')
%% Gradient-based History Matching by finite difference
Jflux=0;
Jwatercut=0;
Jpres=0;

clear Jtotal_Outer_MO Jtotal_MO_OLSPCA update_MO
% global cost function at the reference field
for tt=1:ddd_welldata
    tttt1=(1+(tt-1)*nWI):(nWI+(tt-1)*nWI);
    Jpres=Jpres+(BHPref(:,tt)-BHP(:,tt))'*inv(Rp(tttt1,tttt1))*(BHPref(:,tt)-BHP(:,tt))/2;
    tttt2=(1+(tt-1)*nWP):(nWP+(tt-1)*nWP);  
    Jwatercut=Jwatercut+(Watercutref(:,tt)-Watercut(:,tt))'*inv(Rfw(tttt2,tttt2))*(Watercutref(:,tt)-Watercut(:,tt))/2; 
    Jflux=Jflux+(Fluxref(:,tt)-Flux(:,tt))'*inv(Rq(tttt2,tttt2))*(Fluxref(:,tt)-Flux(:,tt))/2;
end    
Jtotal=Jflux+Jwatercut+Jpres;
Jref_L=Jtotal;
Jref_L

% local cost function at the reference field
erro1=1e-5;  % iteration criteration
erro2=1e-4;
IteNum=100;
sigma1=1;
sigma2=1;

Vandermonde_Chebyshev=3;    % 1 means Vandermonde perturbation, 2 means Chebyshev perturbation
DP_IP=1;    %  1 means dependent perturbation, 2 means independent perturbation
FOM_ROM=2;                 % 1 means FOM, 2 means ROM for cost function calculation
Weighting_Pattern=2;       % 1 means using sin^2+cos^2 function,  2 means using 1/(1+e^x) + e^x/(1+e^x); 3 means 1/(1+x^2) + x^2/(1+x^2)
OPtimization_Algorithm=1;  % 1 means steepset descent, 2 means BFGS, 3 means Gaussian_Newton

NN=0; XXX=[2 5 8];
% Domain_Dx=[3 4 5 2 8];    Domain_Dy=[3 4 5 8 2];  
Domain_Dx=[3  4  5];    Domain_Dy=[3  4  5];   

EnK_num=XXX(2);

Num_K=0;
for z=1:nz
    Num_K=Num_K+size(Layer{z}.fala,2);
end
Perturb_Matrix_MO= (rand(5*nz*EnK_num,Num_K)-0.5*ones(5*nz*EnK_num,Num_K))/100;

Constant_weights=0.25;   % etting constant weightings
for IIJ=1:2
    
    Jtotal_G=[];
    Weightings=[];
    
    AAA=['------ Multi-objective:' num2str(IIJ)  '------'];
    disp(AAA)       

    EnK_num=XXX(2);
    
    % perturbation of global PCA coefficients
    verify_Linearization=zeros(Num_K,1);
    
%     Perturb_Matrix= (rand(5*EnK_num,Num_K)-0.5*ones(5*EnK_num,Num_K))/100;
    num_peturb=size(Perturb_Matrix_MO,1)+1;
    NN=0;

    for ROLM=1:15

        AAA=['------ Outer-loop:' num2str(ROLM)  '------'];
        disp(AAA)  
        
        num_peturb=size(Perturb_Matrix_MO,1)+1;
        % constructing linear model
        K_Layer=[];
        for jj=1:num_peturb

            if jj<num_peturb
               verify2=verify_Linearization+Perturb_Matrix_MO(jj,:)';
            else
               verify2=verify_Linearization;
            end   

            for z=1:nz
                xx=0;
                for zz=1:z
                    xx=xx+size(Layer{zz}.fala,2);
                end
                verify1=verify2(xx-size(Layer{zz}.fala,2)+1:xx);
                K=Layer{z}.mPermeablility+Layer{z}.fala*verify1;
                [K,S]= OPCA_TwoFacies(Layer{z}.fala, verify1, Layer{z}.mPermeablility, gammaF);
                K=K_mud*exp(log(K_sand/K_mud)*K); 
                PORO=0.25 * (K ./200).^0.1;
                ROCK=convertFrom(K,milli*darcy);
                rock{jj}.perm(1+(z-1)*nx*ny:z*nx*ny,1)=ROCK;
                rock{jj}.poro(1+(z-1)*nx*ny:z*nx*ny,1)=PORO;   
            end

        end 

        Gradient_Calculation_finitedifference;
        
        Flux_ROLM=FluxFD(1+(num_peturb-1)*nWP:end,:);
        Watercut_ROLM=WatercutFD(1+(num_peturb-1)*nWP:end,:);
        BHP_ROLM=BHPFD(1+(num_peturb-1)*nWI:end,:);
        
        for JJI=1:size(Domain_Dx,2)
            
            Domain_Numx=Domain_Dx(JJI);
            Domain_Numy=Domain_Dy(JJI); 

            Initilization_DD1;
            parameterization_new1;
            
            Num_K1=0;
            for z=1:nz
                Num_K1=Num_K1+sum(Layer{z}.falaK);
            end
            
            if Num_K< Num_K1
               XX=pinv(TM_GL)*(Perturb_Matrix_MO)';
               Perturb_Matrix1=XX';
            else
                XX= TM_GL'*TM_GL;
                YY= TM_GL'* (Perturb_Matrix_MO)';
                ZZ=inv(XX)*YY;
                Perturb_Matrix1=ZZ';
            end
            
            Gra_Flux1=zeros(nWP,Num_K1,ddd_welldata);
            for j=1:nWP   % for well-level
                for i=1:ddd_welldata
                    B=zeros(Num_K1,1);
                    % extract the local full column rank perturbation matrix
                    [BB]=RBF(WellP_Domain(j),Layer{1}.falaK,Domain_index,Domain_Numx,Domain_Numy,nz);
                    AA=Perturb_Matrix1(:,BB);
                    b=zeros(length(AA(1,:)),1);
                    A=zeros(length(AA(1,:)),length(AA(1,:)));
                    for k=1:num_peturb-1
                        b=b+(FluxFD(j+(k-1)*nWP,i)-FluxFD(j+(num_peturb-1)*nWP,i))*AA(k,:)';
                        A=A+AA(k,:)'*AA(k,:);
                    end 
                    Gra_Flux1(j,BB,i)=pinv(A)*b;
                end
            end

            if size(TM_GL,1)< size(TM_GL,2)
               xx=TM_GL*TM_GL';
               xx=TM_GL'*inv(xx);
            else
                xx=pinv(TM_GL);
            end
            
            for j=1:ddd_welldata
                Multi_Objective{JJI}.Gra_Flux(1:nWP,1:Num_K,j)= Gra_Flux1(:,:,j)*xx;
            end

            Gra_Wat1=zeros(nWP,Num_K1,ddd_welldata);
            for j=1:nWP   % for well-level
                for i=1:ddd_welldata
                    B=zeros(Num_K1,1);
                    % extract the local full column rank perturbation matrix
                    [BB]=RBF(WellP_Domain(j),Layer{1}.falaK,Domain_index,Domain_Numx,Domain_Numy,nz);
                    AA=Perturb_Matrix1(:,BB);
                    b=zeros(length(AA(1,:)),1);
                    A=zeros(length(AA(1,:)),length(AA(1,:)));
                    for k=1:num_peturb-1
                        b=b+(WatercutFD(j+(k-1)*nWP,i)-WatercutFD(j+(num_peturb-1)*nWP,i))*AA(k,:)';
                        A=A+AA(k,:)'*AA(k,:);
                    end 
                    Gra_Wat1(j,BB,i)=pinv(A)*b;
                end
            end

            for j=1:ddd_welldata
                Multi_Objective{JJI}.Gra_Wat(1:nWP,1:Num_K,j)= Gra_Wat1(:,:,j)*xx;
            end

            Gra_BHP1=zeros(nWI,Num_K1,ddd_welldata);
            for j=1:nWI   % for well-level
                for i=1:ddd_welldata
                    B=zeros(Num_K1,1);
                    % extract the local full column rank perturbation matrix
                    [BB]=RBF(WellI_Domain(j),Layer{1}.falaK,Domain_index,Domain_Numx,Domain_Numy,nz);
                    AA=Perturb_Matrix1(:,BB);
                    b=zeros(length(AA(1,:)),1);
                    A=zeros(length(AA(1,:)),length(AA(1,:)));
                    for k=1:num_peturb-1
                        b=b+(BHPFD(j+(k-1)*nWI,i)-BHPFD(j+(num_peturb-1)*nWI,i))*AA(k,:)';
                        A=A+AA(k,:)'*AA(k,:);
                    end 
                    Gra_BHP1(j,BB,i)=pinv(A)*b;
                end
            end 

            for j=1:ddd_welldata
                Multi_Objective{JJI}.Gra_BHP(1:nWI,1:Num_K,j)= Gra_BHP1(:,:,j)*xx;
            end
        
        end
        
        num_peturb=1;
        if ROLM==1
           verify=zeros(Num_K,1);
%            verify=(rand(Num_K,1)-0.5*ones(Num_K,1))*1.5;
           verify_init=verify;
           if length(Domain_Dx)==2
               if IIJ==1 
                   cata=pi/4;
               elseif IIJ==2 && length(Domain_Dx)==2
                   cata=0;
               elseif IIJ==3
                   cata=1;
               else
                   cata=0;
               end
               
           else
               
               if IIJ==1 
                   cata=0;
               elseif IIJ==2
                   cata=1;
               end
                     
           end
           
            xyz=FOM_ROM;
            FOM_ROM=1;

            for z=1:nz
                xx=0;
                for zz=1:z
                    xx=xx+size(Layer{zz}.fala,2);
                end
                verify1=verify(xx-size(Layer{zz}.fala,2)+1:xx);
                K=Layer{z}.mPermeablility+Layer{z}.fala*verify1;

                [K,S]= OPCA_TwoFacies(Layer{z}.fala, verify1, Layer{z}.mPermeablility, gammaF);
                K=K_mud*exp(log(K_sand/K_mud)*K); 
                
                PORO=0.25 * (K ./200).^0.1;
                ROCK=convertFrom(K,milli*darcy);
                rock{1}.perm(1+(z-1)*nx*ny:z*nx*ny,1)=ROCK;
                rock{1}.poro(1+(z-1)*nx*ny:z*nx*ny,1)=PORO;  
            end

            Simulation_FOM_ROM;

            for i=1:num_peturb
                Jflux=0;
                Jwatercut=0;
                Jpres=0;
                tttt1=(1+(i-1)*nWI):(nWI+(i-1)*nWI);
                tttt2=(1+(i-1)*nWP):(length(rSol.wellSol)-nWI+(i-1)*nWP);
                for tt=1:ddd_welldata
                    tttt3=(1+(tt-1)*nWI):(nWI+(tt-1)*nWI);
                    Jpres=Jpres+(BHP(:,tt)-BHPFD(tttt1,tt))'*inv(Rp(tttt3,tttt3))*(BHP(:,tt)-BHPFD(tttt1,tt))/2;
                    tttt4=(1+(tt-1)*nWP):(nWP+(tt-1)*nWP);
                    Jflux=Jflux+(Flux(:,tt)-FluxFD(tttt2,tt))'*inv(Rq(tttt4,tttt4))*(Flux(:,tt)-FluxFD(tttt2,tt))/2;
                    Jwatercut=Jwatercut+(Watercut(:,tt)-WatercutFD(tttt2,tt))'*inv(Rfw(tttt4,tttt4))*(Watercut(:,tt)-WatercutFD(tttt2,tt))/2;         
                end  
                Jtotal=Jflux+Jwatercut+Jpres;
                Jflux=0;
                Jwatercut=0;
                Jpres=0;
            end
            Jtotal_Outer_MO{IIJ}.CF(1,xyz)=Jtotal;
            Jtotal_Outer_MO{IIJ}.Index(1,xyz)=1;
            Jtotal
            
            FOM_ROM=xyz;

        end

        Stop_Condition=0;
        for JJ=1:IteNum  %IteNum
            
            tic 
            if length(Domain_Dx)==2
                
                if IIJ==1
                   weightings=[sin(cata)^2  cos(cata)^2];
                elseif IIJ==2
                   weightings=[1/(1+exp(cata))  exp(cata)/(1+exp(cata))];
                elseif IIJ==3
                   weightings=[1/(1+cata^2)  cata^2/(1+cata^2)];
                else
                   weightings=[Constant_weights  1-Constant_weights];
                end
                
            elseif length(Domain_Dx)>2
                
                if IIJ==1
                    xx=0;
                   for i=1:length(Domain_Dx)
                       xx=xx+exp(cata*(i-1));
                   end
                   for i=1:length(Domain_Dx)
                       weightings(i)= exp(cata*(i-1))/xx;
                   end
                   
                elseif IIJ==2
                    
                   xx=0;
                   for i=1:length(Domain_Dx)
                       xx=xx+cata^((i-1)*2);
                   end
                   for i=1:length(Domain_Dx)
                       weightings(i)= cata^((i-1)*2)/xx;
                   end
                   
                end              
                
            end
            
            for z=1:nz
                xx=0;
                for zz=1:z
                    xx=xx+size(Layer{zz}.fala,2);
                end
                verify1=verify(xx-size(Layer{zz}.fala,2)+1:xx);
                K=Layer{z}.mPermeablility+Layer{z}.fala*verify1;
                [K,S]= OPCA_TwoFacies(Layer{z}.fala, verify1, Layer{z}.mPermeablility, gammaF);
                K=K_mud*exp(log(K_sand/K_mud)*K); 
                PORO=0.25 * (K ./200).^0.1;
                ROCK=convertFrom(K,milli*darcy);
                rock{1}.perm(1+(z-1)*nx*ny:z*nx*ny,1)=ROCK;
                rock{1}.poro(1+(z-1)*nx*ny:z*nx*ny,1)=PORO;   
            end

            if IIJ==1 && JJ==1 && ROLM==1
               K_init=log(convertTo(rock{1}.perm,milli*darcy));
            end
            
            Jtotal11=0;
            Jtotal_weighting=[];
            for JJI=1:size(Domain_Dx,2)
                
                Gra_Flux2=Multi_Objective{JJI}.Gra_Flux;
                Gra_Wat2=Multi_Objective{JJI}.Gra_Wat;
                Gra_BHP2=Multi_Objective{JJI}.Gra_BHP;
                
                Simulation_FOM_ROM_MOT;
 
                % global cost function
                JtotalFDm_G=[];
                for i=1:1
                    Jflux=0;
                    Jwatercut=0;
                    Jpres=0;
                    tttt1=(1+(i-1)*nWI):(nWI+(i-1)*nWI);
                    tttt2=(1+(i-1)*nWP):(length(rSol.wellSol)-nWI+(i-1)*nWP);
                    for tt=1:ddd_welldata
                        tttt3=(1+(tt-1)*nWI):(nWI+(tt-1)*nWI);
                        Jpres=Jpres+(BHP(:,tt)-BHPFD(tttt1,tt))'*inv(Rp(tttt3,tttt3))*(BHP(:,tt)-BHPFD(tttt1,tt))/2;
                        tttt4=(1+(tt-1)*nWP):(nWP+(tt-1)*nWP);
                        Jflux=Jflux+(Flux(:,tt)-FluxFD(tttt2,tt))'*inv(Rq(tttt4,tttt4))*(Flux(:,tt)-FluxFD(tttt2,tt))/2;
                        Jwatercut=Jwatercut+(Watercut(:,tt)-WatercutFD(tttt2,tt))'*inv(Rfw(tttt4,tttt4))*(Watercut(:,tt)-WatercutFD(tttt2,tt))/2;         
                    end  
                    Jtotal=Jflux+Jwatercut+Jpres;
                    Jflux=0;
                    Jwatercut=0;
                    Jpres=0;
                end
                Jtotal11=Jtotal11+Jtotal*weightings(JJI);
                Jtotal_weighting=[Jtotal_weighting Jtotal];
                
            end
            Jtotal11
            Jtotal_G=[Jtotal_G Jtotal11];
            
            if JJ>1   

                sigma1=abs(Jtotal_G(end)-Jtotal_G(end-1))/max([1 Jtotal_G(end)]);
                sigma2=norm(verify-verifym,2)/norm(verify,2);

                if sigma1<erro1 
                   Stop_Condition=1;
                elseif sigma2<erro2
                   Stop_Condition=1;
                elseif (JJ-1)>=IteNum
                   Stop_Condition=1;
                end

            end

            if Stop_Condition==1

               break
               
            end

            EE=verify;
            for JJI=1:size(Domain_Dx,2)
                
                Gra_Flux2=Multi_Objective{JJI}.Gra_Flux;
                Gra_Wat2=Multi_Objective{JJI}.Gra_Wat;
                Gra_BHP2=Multi_Objective{JJI}.Gra_BHP;
                
                Simulation_FOM_ROM_MOT;
                
                Gradient_Computation_MO;
                
            end
            
            if length(Domain_Dx)==2
                
                if IIJ==1
                   EE=[EE ; sin(2*cata)*(Jtotal_weighting(1)-Jtotal_weighting(2))];
                elseif IIJ==2
                   EE=[EE ; exp(cata)/(1+exp(cata))^2*(Jtotal_weighting(2)-Jtotal_weighting(1))];
                elseif IIJ==3
                   EE=[EE ; 2*cata/(1+cata^2)^2*(Jtotal_weighting(2)-Jtotal_weighting(1))];
                else
                   EE=[EE ; 0];
                end
                
            elseif length(Domain_Dx)>2
                
                if IIJ==1
                   xx=0;
                   for i=1:length(Domain_Dx)
                       xx=xx+exp(cata*(i-1));
                   end
                   
                   yy=0;
                   for i=1:length(Domain_Dx)
                       yy=yy+(i-1)*exp(cata*(i-1));
                   end
                   
                   zz=0;
                   for i=1:length(Domain_Dx)
                       zz=zz+((i-1)*exp(cata*(i-1))*xx-exp(cata*(i-1))*yy)/xx/xx*Jtotal_weighting(i);
                   end
                   EE=[EE ; zz];
                   
                elseif IIJ==2
                    
                   xx=0;
                   for i=1:length(Domain_Dx)
                       xx=xx+cata^((i-1)*2);
                   end
                   
                   yy=0;
                   for i=1:length(Domain_Dx)
                       yy=yy+2*(i-1)*cata^((i-1)*2-1);
                   end
                   
                   zz=0;
                   for i=1:length(Domain_Dx)
                       zz=zz+(2*(i-1)*cata^((i-1)*2-1)*xx-cata^((i-1)*2)*yy)/xx/xx*Jtotal_weighting(i);
                   end
                   EE=[EE ; zz];
                   
                end
                
            end
            
            EE=EE/max(abs(EE));
            
            
            if OPtimization_Algorithm<3
           
                if JJ==1 
                   ll=0.2;
                end

               if  JJ>1 && Jtotal_G(end)>Jtotal_G(end-1)
                    ll=ll/2;
                    verify=verifym-ll*FF(1:end-1);
                    cata=catam-ll*FF(end);
                else 

                     verifym=verify;
                     catam=cata;
                     verify=verifym-ll*EE(1:end-1);  
                     cata=catam-ll*EE(end);
               end
                
            elseif OPtimization_Algorithm == 3
                
                   verifym=verify;
                   verify=verifym-EE;  
                
            end

            FF=EE;
            
           cata
           Weightings=[Weightings cata];

            AAA=['Inner-loop Iteration Step:' num2str(JJ)];
            disp(AAA)

            toc
            
        end
        Stop_Condition=0;

        for z=1:nz
            xx=0;
            for zz=1:z
                xx=xx+size(Layer{zz}.fala,2);
            end
            verify1=verify(xx-size(Layer{zz}.fala,2)+1:xx);
            K_upt=Layer{z}.mPermeablility+Layer{z}.fala*verify1;

            [K_upt,S]= OPCA_TwoFacies(Layer{z}.fala, verify1, Layer{z}.mPermeablility, gammaF);
            K=K_mud*exp(log(K_sand/K_mud)*K_upt); 

            PORO=0.25 * (K ./200).^0.1;
            ROCK=convertFrom(K,milli*darcy);
            rock{1}.perm(1+(z-1)*nx*ny:z*nx*ny,1)=ROCK;
            rock{1}.poro(1+(z-1)*nx*ny:z*nx*ny,1)=PORO; 
        end       
        update_MO{IIJ}.Perm(1:length(K_upt),ROLM)= K_upt;
        
        verify_Linearization=verify;
        
        
        xyz=FOM_ROM;
        FOM_ROM=1;

        Simulation_FOM_ROM;

        for i=1:num_peturb
            Jflux=0;
            Jwatercut=0;
            Jpres=0;
            tttt1=(1+(i-1)*nWI):(nWI+(i-1)*nWI);
            tttt2=(1+(i-1)*nWP):(length(rSol.wellSol)-nWI+(i-1)*nWP);
            for tt=1:ddd_welldata
                tttt3=(1+(tt-1)*nWI):(nWI+(tt-1)*nWI);
                Jpres=Jpres+(BHP(:,tt)-BHPFD(tttt1,tt))'*inv(Rp(tttt3,tttt3))*(BHP(:,tt)-BHPFD(tttt1,tt))/2;
                tttt4=(1+(tt-1)*nWP):(nWP+(tt-1)*nWP);
                Jflux=Jflux+(Flux(:,tt)-FluxFD(tttt2,tt))'*inv(Rq(tttt4,tttt4))*(Flux(:,tt)-FluxFD(tttt2,tt))/2;
                Jwatercut=Jwatercut+(Watercut(:,tt)-WatercutFD(tttt2,tt))'*inv(Rfw(tttt4,tttt4))*(Watercut(:,tt)-WatercutFD(tttt2,tt))/2;         
            end  
            Jtotal=Jflux+Jwatercut+Jpres;
            Jflux=0;
            Jwatercut=0;
            Jpres=0;
        end
        Jtotal_Outer_MO{IIJ}.CF(ROLM+1,xyz)=Jtotal;
        Jtotal_Outer_MO{IIJ}.Index(ROLM+1,xyz)=length(Jtotal_G);
        Jtotal


        FOM_ROM=xyz;


    end
    Jtotal_MO_OLSPCA{IIJ}.G(1:length(Jtotal_G),FOM_ROM)=Jtotal_G';
    
    
    if length(Domain_Dx)==2

        if IIJ==1
           
            for i=1:length(Weightings)
                Jtotal_MO_OLSPCA{IIJ}.weights(i,1)=(sin(Weightings(i)))^2;
                Jtotal_MO_OLSPCA{IIJ}.weights(i,2)=(cos(Weightings(i)))^2;
            end
        elseif IIJ==2
            
            for i=1:length(Weightings)
                Jtotal_MO_OLSPCA{IIJ}.weights(i,1)=1/(1+exp(Weightings(i)));
                Jtotal_MO_OLSPCA{IIJ}.weights(i,2)=exp(Weightings(i))/(1+exp(Weightings(i)));
            end 

        elseif IIJ==3
            
            for i=1:length(Weightings)
                Jtotal_MO_OLSPCA{IIJ}.weights(i,1)=1/(1+(Weightings(i))^2);
                Jtotal_MO_OLSPCA{IIJ}.weights(i,2)=(Weightings(i))^2/(1+(Weightings(i))^2);
            end 

        else
            
            for i=1:length(Weightings)
                Jtotal_MO_OLSPCA{IIJ}.weights(i,1)=Constant_weights;
                Jtotal_MO_OLSPCA{IIJ}.weights(i,2)=1-Constant_weights;
            end 

        end

    elseif length(Domain_Dx)>2

        if IIJ==1
            
            for j=1:length(Weightings)
                
                xx=0;
               for i=1:length(Domain_Dx)
                   xx=xx+exp(Weightings(j)*(i-1));
               end
               for i=1:length(Domain_Dx)
                   Jtotal_MO_OLSPCA{IIJ}.weights(j,i)= exp(Weightings(j)*(i-1))/xx;
               end
               
            end

        elseif IIJ==2

            for j=1:length(Weightings)
                
               xx=0;
               for i=1:length(Domain_Dx)
                   xx=xx+(Weightings(j))^((i-1)*2);
               end
               for i=1:length(Domain_Dx)
                   Jtotal_MO_OLSPCA{IIJ}.weights(j,i)= (Weightings(j))^((i-1)*2)/xx;
               end
           
            end

        end              

    end
    

    figure(1)
    subplot(1,2,1)
    semilogy(Jtotal_G); hold on
    if IIJ==1
        semilogy(1:IteNum*ROLM,ones(IteNum*ROLM,1)*(2*nWP+nWI)*ddd*5/2,'k--','LineWidth',2);
        semilogy(1:IteNum*ROLM,ones(IteNum*ROLM,1)*Jref_L,'r--','LineWidth',2);
    end
    
    subplot(1,2,2)
    semilogy(Jtotal_Outer_MO{IIJ}.CF(:,xyz),'-o','LineWidth',1); hold on
    if IIJ==1
        semilogy(1:ROLM+1,ones(1+ROLM,1)*(2*nWP+nWI)*ddd*5/2,'k--','LineWidth',2); 
        semilogy(1:1+ROLM,ones(1+ROLM,1)*Jref_L,'r--','LineWidth',2);
    end
    
    figure(2)
    for j=1:length(Domain_Dx)
        plot(Jtotal_MO_OLSPCA{IIJ}.weights(:,j)); hold on
    end

        
end

