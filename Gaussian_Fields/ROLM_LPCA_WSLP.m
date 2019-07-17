%% Gradient-based History Matching by finite difference
Jseismic=0;
Jflux=0;
Jwatercut=0;
Jpres=0;
clear Multi_Objective
clear Jtotal_Outer_MO_SeismicData_ProductionData Jtotal_MO_OLSPCA_SeismicData_ProductionData update_MO_SeismicData_ProductionData
% global cost function at the reference field
for tt=1:ddd_welldata
    tttt1=(1+(tt-1)*nWI):(nWI+(tt-1)*nWI);
    Jpres=Jpres+(BHPref(:,tt)-BHP(:,tt))'*inv(Rp(tttt1,tttt1))*(BHPref(:,tt)-BHP(:,tt))/2;
    tttt2=(1+(tt-1)*nWP):(nWP+(tt-1)*nWP);  
    Jwatercut=Jwatercut+(Watercutref(:,tt)-Watercut(:,tt))'*inv(Rfw(tttt2,tttt2))*(Watercutref(:,tt)-Watercut(:,tt))/2; 
    Jflux=Jflux+(Fluxref(:,tt)-Flux(:,tt))'*inv(Rq(tttt2,tttt2))*(Fluxref(:,tt)-Flux(:,tt))/2;
end    
Jtotal=Jflux+Jwatercut+Jpres;

Jseismic=0;
for tt=1:ddd_seismicdata
    tttt1=(1+(tt-1)*nx*ny*nz):(tt*nx*ny*nz);
    Jseismic=Jseismic+(Satref(:,tt)-Sat(:,tt))'*inv(Rs(tttt1,tttt1))*(Satref(:,tt)-Sat(:,tt))/2;
end  
Jref_L=Jtotal+Jseismic;
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


Num_kr=6;     % updating 6 relative permeability parameters
Both_K_Kr=2;  % 1 only updating log-permeability, 2 simultaneously updating 6 relative permeability parameters

Domain_Dx=[3  4  5];    Domain_Dy=[4 5 6]; 

Domain_Numx=min(Domain_Dx);
Domain_Numy=min(Domain_Dy); 

Initilization_DD1;

Num_LPCA=(ceil(size(Layer{nz}.fala,2)/Domain_Numx/Domain_Numy))*ones(Domain_Numx*Domain_Numy,1);

parameterization_new1;

Num_K=0;
for z=1:nz
    Num_K=Num_K+size(Layer{z}.fala,2);
end
Perturb_Matrix_MO= (rand(5*nz*Num_LPCA(1),Num_K)-0.5*ones(5*nz*Num_LPCA(1),Num_K))/100;

Perturb_Matrix=Perturb_Matrix_MO;

Constant_weights=0.25;   % etting constant weightings
for IIJ=1:1
    
    Jtotal_G=[];
    Weightings=[];
    Weightings_Out=[];
      
    AAA=['------ Multi-objective:' num2str(IIJ)  '------'];
    disp(AAA)       

    % perturbation of global PCA coefficients
    verify_Linearization = zeros(Num_K+Num_kr,1); 
    verify_kr_initial=rand(Num_kr,1);
    verify_Linearization(Num_K+1:Num_K+Num_kr,1)=verify_kr_initial;

    Perturb_Matrix_kr=delta*eye(Num_kr,Num_kr);
    
    NN=0;

    for ROLM=1:10

        AAA=['------ Outer-loop:' num2str(ROLM)  '------'];
        disp(AAA)  
        
        num_peturb=size(Perturb_Matrix_MO,1)+size(Perturb_Matrix_kr,1)+1;
        % constructing linear model
        K_Layer=[];
        for jj=1:num_peturb

            if jj<size(Perturb_Matrix_MO,1)+1
               verify2=verify_Linearization(1:Num_K)+Perturb_Matrix_MO(jj,:)';
            else
               verify2=verify_Linearization(1:Num_K);
            end   

            for z=1:nz
                xx=0;
                for zz=1:z
                    xx=xx+size(Layer{z}.fala,2);
                end
                verify1=verify2(xx-size(Layer{z}.fala,2)+1:xx);
                K=Layer{z}.mPermeablility+Layer{z}.fala*verify1;
                PORO=0.25 * (exp(K) ./200).^0.1;
                ROCK=convertFrom(exp(K),milli*darcy);
                rock{jj}.perm(1+(z-1)*nx*ny:z*nx*ny,1)=ROCK;
                rock{jj}.poro(1+(z-1)*nx*ny:z*nx*ny,1)=PORO;    
            end

        end 

        Gradient_Calculation_finitedifference;
        
        Flux_ROLM=FluxFD(1+(num_peturb-1)*nWP:end,:);
        Watercut_ROLM=WatercutFD(1+(num_peturb-1)*nWP:end,:);
        BHP_ROLM=BHPFD(1+(num_peturb-1)*nWI:end,:);
        Seismic_ROLM=SatFD(1+(num_peturb-1)*nx*ny*nz:end,:);
        
        for JJI=1:size(Domain_Dx,2)
            
            Domain_Numx=Domain_Dx(JJI);
            Domain_Numy=Domain_Dy(JJI); 

            Initilization_DD1;
            
            Num_LPCA=(ceil(size(Layer{nz}.fala,2)/Domain_Numx/Domain_Numy))*ones(Domain_Numx*Domain_Numy,1);
            
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
            
            Gra_Seismic1=zeros(nx*ny*nz,Num_K1,ddd_seismicdata);
            for j=1:Domain_Numx*Domain_Numy   % for well-level

                Position1=[];
                for ij=Domain_index(j,5):Domain_index(j,6)
                    for kk=Domain_index(j,3):Domain_index(j,4)
                        Position1=[Position1 kk+(ij-1)*nx];
                    end
                end
                Position=[];
                for z=1:nz
                    Position=[Position Position1+ones(size(Position1))*(z-1)*nx*ny];
                end

                for i=1:ddd_seismicdata
                    % extract the local full column rank perturbation matrix
                    [BB]=RBF(j,Layer{1}.falaK,Domain_index,Domain_Numx,Domain_Numy,nz);
                    AA=Perturb_Matrix1(:,BB);
                    for z=1:length(Position)
                        b=zeros(length(AA(1,:)),1);
                        A=zeros(length(AA(1,:)),length(AA(1,:)));
                        for k=1:size(Perturb_Matrix1,1)
                            b=b+(SatFD(Position(z)+(k-1)*nx*ny*nz,i)-SatFD(Position(z)+(num_peturb-1)*nx*ny*nz,i))*AA(k,:)';
                            A=A+AA(k,:)'*AA(k,:);
                        end
                        Gra_Seismic1(Position(z),BB,i)=pinv(A)*b;
                    end
                end
            end
            
            if size(TM_GL,1)< size(TM_GL,2)
               xx=TM_GL*TM_GL';
               xx=TM_GL'*inv(xx);
            else
                xx=pinv(TM_GL);
            end
            
            for j=1:ddd_seismicdata
                Multi_Objective{JJI}.Gra_Seismic(1:nx*ny*nz,1:Num_K,j)= Gra_Seismic1(:,:,j)*xx;
            end
            
            
            for j=1:nx*ny*nz   % for well-level
                for i=1:ddd_seismicdata
                    for k=size(Perturb_Matrix1,1)+1:num_peturb-1
                        Multi_Objective{JJI}.Gra_Seismic(j,k-size(Perturb_Matrix1,1)+Num_K,i)= (SatFD(j+(k-1)*nx*ny*nz,i)-SatFD(j+(num_peturb-1)*nx*ny*nz,i))/0.01; 
                    end  
                end
            end
            
            Gra_Flux1=zeros(nWP,Num_K1,ddd_welldata);
            for j=1:nWP   % for well-level
                for i=1:ddd_welldata
                    % extract the local full column rank perturbation matrix
                    if Vandermonde_Chebyshev==3
                       [BB]=RBF(WellP_Domain(j),Layer{1}.falaK,Domain_index,Domain_Numx,Domain_Numy,nz);
                    else
                        BB=[];
                        for z=1:nz
                            tt1=sum(Layer{z}.falaK(1:WellP_Domain(j)))-Layer{z}.falaK(WellP_Domain(j))+1:sum(Layer{z}.falaK(1:WellP_Domain(j)));
                            BB=[BB tt1];
                        end  
                    end

                    AA=Perturb_Matrix1(:,BB);
                    b=zeros(length(AA(1,:)),1);
                    A=zeros(length(AA(1,:)),length(AA(1,:)));
                    for k=1:size(Perturb_Matrix1,1)
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
            
            
            for j=1:nWP   % for well-level
                for i=1:ddd_welldata
                    for k=size(Perturb_Matrix1,1)+1:num_peturb-1
                        Multi_Objective{JJI}.Gra_Flux(j,k-size(Perturb_Matrix1,1)+Num_K,i)= (FluxFD(j+(k-1)*nWP,i)-FluxFD(j+(num_peturb-1)*nWP,i))/0.01; 
                    end  
                end
            end
            

            Gra_Wat1=zeros(nWP,Num_K1,ddd_welldata);
            for j=1:nWP   % for well-level
                for i=1:ddd_welldata
                    % extract the local full column rank perturbation matrix
                    if Vandermonde_Chebyshev==3
                       [BB]=RBF(WellP_Domain(j),Layer{1}.falaK,Domain_index,Domain_Numx,Domain_Numy,nz);
                    else
                        BB=[];
                        for z=1:nz
                            tt1=sum(Layer{z}.falaK(1:WellP_Domain(j)))-Layer{z}.falaK(WellP_Domain(j))+1:sum(Layer{z}.falaK(1:WellP_Domain(j)));
                            BB=[BB tt1];
                        end  
                    end
                    AA=Perturb_Matrix1(:,BB);
                    b=zeros(length(AA(1,:)),1);
                    A=zeros(length(AA(1,:)),length(AA(1,:)));
                    for k=1:size(Perturb_Matrix1,1)
                        b=b+(WatercutFD(j+(k-1)*nWP,i)-WatercutFD(j+(num_peturb-1)*nWP,i))*AA(k,:)';
                        A=A+AA(k,:)'*AA(k,:);
                    end 
                    Gra_Wat1(j,BB,i)=pinv(A)*b;
                end
            end
            
            if size(TM_GL,1)< size(TM_GL,2)
               xx=TM_GL*TM_GL';
               xx=TM_GL'*inv(xx);
            else
                xx=pinv(TM_GL);
            end
            
            for j=1:ddd_welldata
                Multi_Objective{JJI}.Gra_Wat(1:nWP,1:Num_K,j)= Gra_Wat1(:,:,j)*xx;
            end
            
            
            for j=1:nWP   % for well-level
                for i=1:ddd_welldata
                    for k=size(Perturb_Matrix1,1)+1:num_peturb-1
                        Multi_Objective{JJI}.Gra_Wat(j,k-size(Perturb_Matrix1,1)+Num_K,i)= (WatercutFD(j+(k-1)*nWP,i)-WatercutFD(j+(num_peturb-1)*nWP,i))/0.01; 
                    end  
                end
            end

            Gra_BHP1=zeros(nWI,Num_K1,ddd_welldata);
            for j=1:nWI   % for well-level
                for i=1:ddd_welldata
                    % extract the local full column rank perturbation matrix
                    if Vandermonde_Chebyshev==3
                       [BB]=RBF(WellI_Domain(j),Layer{1}.falaK,Domain_index,Domain_Numx,Domain_Numy,nz);
                    else
                        BB=[];
                        for z=1:nz
                            tt1=sum(Layer{z}.falaK(1:WellI_Domain(j)))-Layer{z}.falaK(WellI_Domain(j))+1:sum(Layer{z}.falaK(1:WellI_Domain(j)));
                            BB=[BB tt1];
                        end  
                    end
                    AA=Perturb_Matrix1(:,BB);
                    b=zeros(length(AA(1,:)),1);
                    A=zeros(length(AA(1,:)),length(AA(1,:)));
                    for k=1:size(Perturb_Matrix1,1)
                        b=b+(BHPFD(j+(k-1)*nWI,i)-BHPFD(j+(num_peturb-1)*nWI,i))*AA(k,:)';
                        A=A+AA(k,:)'*AA(k,:);
                    end 
                    Gra_BHP1(j,BB,i)=pinv(A)*b;
                end
            end
            
            if size(TM_GL,1)< size(TM_GL,2)
               xx=TM_GL*TM_GL';
               xx=TM_GL'*inv(xx);
            else
                xx=pinv(TM_GL);
            end
            
            for j=1:ddd_welldata
                Multi_Objective{JJI}.Gra_BHP(1:nWI,1:Num_K,j)= Gra_BHP1(:,:,j)*xx;
            end
            
            
            for j=1:nWI   % for well-level
                for i=1:ddd_welldata
                    for k=size(Perturb_Matrix1,1)+1:num_peturb-1
                        Multi_Objective{JJI}.Gra_BHP(j,k-size(Perturb_Matrix1,1)+Num_K,i)= (BHPFD(j+(k-1)*nWI,i)-BHPFD(j+(num_peturb-1)*nWI,i))/0.01; 
                    end  
                end
            end

        end
        
        num_peturb=1;
        if ROLM==1
            
           verify=verify_Linearization;
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
                PORO=0.25 * (exp(K) ./200).^0.1;
                ROCK=convertFrom(exp(K),milli*darcy);
                rock{1}.perm(1+(z-1)*nx*ny:z*nx*ny,1)=ROCK;
                rock{1}.poro(1+(z-1)*nx*ny:z*nx*ny,1)=PORO;    
            end

            Simulation_FOM_ROM_Seismic_ProductionData;

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
            
            Jseismic=0;
            for tt=1:ddd_seismicdata
                tttt1=(1+(tt-1)*nx*ny*nz):(tt*nx*ny*nz);
                Jseismic=Jseismic+(SatFD(:,tt)-Sat(:,tt))'*inv(Rs(tttt1,tttt1))*(SatFD(:,tt)-Sat(:,tt))/2;
            end 
            
            Jtotal_Outer_MO_SeismicData_ProductionData{IIJ}.CF(1,xyz)=Jtotal+Jseismic;
            Jtotal_Outer_MO_SeismicData_ProductionData{IIJ}.Index(1,xyz)=1;
            Jtotal+Jseismic
            
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
                PORO=0.25 * (exp(K) ./200).^0.1;
                ROCK=convertFrom(exp(K),milli*darcy);
                rock{1}.perm(1+(z-1)*nx*ny:z*nx*ny,1)=ROCK;
                rock{1}.poro(1+(z-1)*nx*ny:z*nx*ny,1)=PORO;    
            end

            if IIJ==1 && JJ==1 && ROLM==1
               K_init=log(convertTo(rock{1}.perm,milli*darcy));
            end
            
            Jtotal11=0;
            Jtotal_weighting=[];
            for JJI=1:size(Domain_Dx,2)
                
                Gra_Seismic2=Multi_Objective{JJI}.Gra_Seismic;
                Gra_BHP2=Multi_Objective{JJI}.Gra_BHP;
                Gra_Wat2=Multi_Objective{JJI}.Gra_Wat;
                Gra_Flux2=Multi_Objective{JJI}.Gra_Flux;
                
                Simulation_FOM_ROM_MOT_Seismic_ProductionData;
 
                % global cost function
                Jseismic=0;
                for tt=1:ddd_seismicdata
                    tttt1=(1+(tt-1)*nx*ny*nz):(tt*nx*ny*nz);
                    Jseismic=Jseismic+(SatFD(:,tt)-Sat(:,tt))'*inv(Rs(tttt1,tttt1))*(SatFD(:,tt)-Sat(:,tt))/2;
                end
                
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
                
                Jtotal11=Jtotal11+(Jseismic+Jtotal)*weightings(JJI);
                Jtotal_weighting=[Jtotal_weighting (Jseismic+Jtotal)];
                
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
                
                Gra_Seismic2=Multi_Objective{JJI}.Gra_Seismic;
                Gra_BHP2=Multi_Objective{JJI}.Gra_BHP;
                Gra_Wat2=Multi_Objective{JJI}.Gra_Wat;
                Gra_Flux2=Multi_Objective{JJI}.Gra_Flux;
                
                Simulation_FOM_ROM_MOT_Seismic_ProductionData;
                
                Gradient_Computation_MO_Seismic_ProductionData;
                
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
        
        Weightings_Out=[Weightings_Out cata];
        
        Stop_Condition=0;

        for z=1:nz
            xx=0;
            for zz=1:z
                xx=xx+size(Layer{zz}.fala,2);
            end
            verify1=verify(xx-size(Layer{zz}.fala,2)+1:xx);
            K_upt=Layer{z}.mPermeablility+Layer{z}.fala*verify1;
            PORO=0.25 * (exp(K_upt) ./200).^0.1;
            ROCK=convertFrom(exp(K_upt),milli*darcy);
            rock{1}.perm(1+(z-1)*nx*ny:z*nx*ny,1)=ROCK;
            rock{1}.poro(1+(z-1)*nx*ny:z*nx*ny,1)=PORO;    
        end       
        update_MO_SeismicData_ProductionData{IIJ}.Perm(1:length(K_upt),ROLM)= K_upt;
        
        verify_Linearization=verify;
        
        
        xyz=FOM_ROM;
        FOM_ROM=1;

        Simulation_FOM_ROM_Seismic_ProductionData;

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

        Jseismic=0;
        for tt=1:ddd_seismicdata
            tttt1=(1+(tt-1)*nx*ny*nz):(tt*nx*ny*nz);
            Jseismic=Jseismic+(SatFD(:,tt)-Sat(:,tt))'*inv(Rs(tttt1,tttt1))*(SatFD(:,tt)-Sat(:,tt))/2;
        end 

        Jtotal_Outer_MO_SeismicData_ProductionData{IIJ}.CF(ROLM+1,xyz)=Jtotal+Jseismic;
        Jtotal_Outer_MO_SeismicData_ProductionData{IIJ}.Index(ROLM+1,xyz)=length(Jtotal_G);
        Jtotal+Jseismic

        FOM_ROM=xyz;
         
    end
    Jtotal_MO_OLSPCA_SeismicData_ProductionData{IIJ}.G(1:length(Jtotal_G),FOM_ROM)=Jtotal_G';

    if length(Domain_Dx)==2

        if IIJ==1
           
            for i=1:length(Weightings)
                Jtotal_MO_OLSPCA_SeismicData_ProductionData{IIJ}.weights(i,1)=sin(Weightings(i));
                Jtotal_MO_OLSPCA_SeismicData_ProductionData{IIJ}.weights(i,2)=cos(Weightings(i));
            end
            
             for i=1:length(Weightings_Out)
                Jtotal_Outer_MO_SeismicData_ProductionData{IIJ}.weights(i,1)=sin(Weightings_Out(i));
                Jtotal_Outer_MO_SeismicData_ProductionData{IIJ}.weights(i,2)=cos(Weightings_Out(i));
            end
            
        elseif IIJ==2
            
            for i=1:length(Weightings)
                Jtotal_MO_OLSPCA_SeismicData_ProductionData{IIJ}.weights(i,1)=1/(1+exp(Weightings(i)));
                Jtotal_MO_OLSPCA_SeismicData_ProductionData{IIJ}.weights(i,2)=exp(Weightings(i))/(1+exp(Weightings(i)));
            end 
            
            for i=1:length(Weightings_Out)
                Jtotal_Outer_MO_SeismicData_ProductionData{IIJ}.weights(i,1)=1/(1+exp(Weightings_Out(i)));
                Jtotal_Outer_MO_SeismicData_ProductionData{IIJ}.weights(i,2)=exp(Weightings_Out(i))/(1+exp(Weightings_Out(i)));
            end

        elseif IIJ==3
            
            for i=1:length(Weightings)
                Jtotal_MO_OLSPCA_SeismicData_ProductionData{IIJ}.weights(i,1)=1/(1+(Weightings(i))^2);
                Jtotal_MO_OLSPCA_SeismicData_ProductionData{IIJ}.weights(i,2)=(Weightings(i))^2/(1+(Weightings(i))^2);
            end 
            
            for i=1:length(Weightings_Out)
                Jtotal_Outer_MO_SeismicData_ProductionData{IIJ}.weights(i,1)=1/(1+(Weightings_Out(i))^2);
                Jtotal_Outer_MO_SeismicData_ProductionData{IIJ}.weights(i,2)=(Weightings_Out(i))^2/(1+(Weightings_Out(i))^2);
            end

        else
            
            for i=1:length(Weightings)
                Jtotal_MO_OLSPCA_SeismicData_ProductionData{IIJ}.weights(i,1)=Constant_weights;
                Jtotal_MO_OLSPCA_SeismicData_ProductionData{IIJ}.weights(i,2)=1-Constant_weights;
            end 
            
            for i=1:length(Weightings_Out)
                Jtotal_Outer_MO_SeismicData_ProductionData{IIJ}.weights(i,1)=Constant_weights;
                Jtotal_Outer_MO_SeismicData_ProductionData{IIJ}.weights(i,2)=1-Constant_weights;
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
                   Jtotal_MO_OLSPCA_SeismicData_ProductionData{IIJ}.weights(j,i)= exp(Weightings(j)*(i-1))/xx;
               end
               
            end
            
            for j=1:length(Weightings_Out)
                
                xx=0;
               for i=1:length(Domain_Dx)
                   xx=xx+exp(Weightings_Out(j)*(i-1));
               end
               for i=1:length(Domain_Dx)
                   Jtotal_Outer_MO_SeismicData_ProductionData{IIJ}.weights(j,i)= exp(Weightings_Out(j)*(i-1))/xx;
               end
               
            end
            

        elseif IIJ==2

            for j=1:length(Weightings)
                
               xx=0;
               for i=1:length(Domain_Dx)
                   xx=xx+(Weightings(j))^((i-1)*2);
               end
               for i=1:length(Domain_Dx)
                   Jtotal_MO_OLSPCA_SeismicData_ProductionData{IIJ}.weights(j,i)= (Weightings(j))^((i-1)*2)/xx;
               end
           
            end
            
            for j=1:length(Weightings_Out)
                
               xx=0;
               for i=1:length(Domain_Dx)
                   xx=xx+(Weightings_Out(j))^((i-1)*2);
               end
               for i=1:length(Domain_Dx)
                   Jtotal_Outer_MO_SeismicData_ProductionData{IIJ}.weights(j,i)= (Weightings_Out(j))^((i-1)*2)/xx;
               end
           
            end

        end              

    end    

    figure(1)
    subplot(1,2,1)
    semilogy(Jtotal_G); hold on
    if IIJ==1
        semilogy(1:1000,ones(1000,1)*size(SatFD,1)*size(SatFD,2)*5,'k--','LineWidth',2); 
        semilogy(1:1000,ones(1000,1)*Jref_L,'r--','LineWidth',2);
    end
    
    subplot(1,2,2)
    semilogy(Jtotal_Outer_MO_SeismicData_ProductionData{IIJ}.CF(:,xyz),'-o','LineWidth',1); hold on
    if IIJ==1
        semilogy(1:ROLM+1,ones(1+ROLM,1)*size(SatFD,1)*size(SatFD,2)*5,'k--','LineWidth',2); 
        semilogy(1:1+ROLM,ones(1+ROLM,1)*Jref_L,'r--','LineWidth',2);
    end
    
    figure(2)
    subplot(1,2,1)
    for j=1:length(Domain_Dx)
        plot(Jtotal_MO_OLSPCA_SeismicData_ProductionData{IIJ}.weights(:,j)); hold on
    end
    
    subplot(1,2,2)
    for j=1:length(Domain_Dx)
        plot(Jtotal_Outer_MO_SeismicData_ProductionData{IIJ}.weights(:,j),'-o'); hold on
    end 
    
    
        
end


figure
fluid = initCoreyFluid('mu' , [ 0.4, 2]*centi*poise     , ...
                       'rho', [1014, 859]*kilogram/meter^3, ...
                       'n'  , [ n_mean, n_mean]                 , ...
                       'sr' , [ S_mean, S_mean]                 , ...
                       'kwm', [ kr_mean,  kr_mean]); 
s = linspace(0, 1, 1001).'; kr = fluid.relperm(s);
plot(s, kr(:,1),'b'); hold on
plot(s, kr(:,2),'b'); hold on


if Num_kr==6 && Both_K_Kr==2

    verify_kr=verify_kr_initial; 

    verify_kr=(exp(verify_kr).*[n_max; n_max; S_max; S_max; kr_max; kr_max]+[n_min; n_min; S_min; S_min; kr_min; kr_min])./(ones(Num_kr,1)+exp(verify_kr));    

    fluid = initCoreyFluid('mu' , [ 0.4, 2]*centi*poise     , ...
                           'rho', [1014, 859]*kilogram/meter^3, ...
                           'n'  , [   verify_kr(1),   verify_kr(2)]                 , ...
                           'sr' , [ verify_kr(3), verify_kr(4)]                 , ...
                           'kwm', [  verify_kr(5),   verify_kr(6)]);

elseif Num_kr==3 && Both_K_Kr==2

    verify_kr=verify_kr_initial; 
    verify_kr=(exp(verify_kr).*[n_max; n_max; kr_max]+[n_min; n_min; kr_min])./(ones(Num_kr,1)+exp(verify_kr));

    fluid = initCoreyFluid('mu' , [ 0.4, 2]*centi*poise     , ...
                           'rho', [1014, 859]*kilogram/meter^3, ...
                           'n'  , [   verify_kr(1),   verify_kr(2)]                 , ...
                           'sr' , [ S_mean, S_mean]                 , ...
                           'kwm', [   verify_kr(3),   kr_mean]);     

end
s = linspace(0, 1, 1001).'; kr = fluid.relperm(s);
plot(s, kr(:,1),'g'); hold on
plot(s, kr(:,2),'g'); hold on


if Num_kr==6 && Both_K_Kr==2

    verify_kr=verify(Num_K+1:Num_K+Num_kr,1); 

    verify_kr=(exp(verify_kr).*[n_max; n_max; S_max; S_max; kr_max; kr_max]+[n_min; n_min; S_min; S_min; kr_min; kr_min])./(ones(Num_kr,1)+exp(verify_kr));    

    fluid = initCoreyFluid('mu' , [ 0.4, 2]*centi*poise     , ...
                           'rho', [1014, 859]*kilogram/meter^3, ...
                           'n'  , [   verify_kr(1),   verify_kr(2)]                 , ...
                           'sr' , [ verify_kr(3), verify_kr(4)]                 , ...
                           'kwm', [  verify_kr(5),   verify_kr(6)]);

elseif Num_kr==3 && Both_K_Kr==2

    verify_kr=verify(Num_K+1:Num_K+Num_kr,1); 
    verify_kr=(exp(verify_kr).*[n_max; n_max; kr_max]+[n_min; n_min; kr_min])./(ones(Num_kr,1)+exp(verify_kr));

    fluid = initCoreyFluid('mu' , [ 0.4, 2]*centi*poise     , ...
                           'rho', [1014, 859]*kilogram/meter^3, ...
                           'n'  , [   verify_kr(1),   verify_kr(2)]                 , ...
                           'sr' , [ S_mean, S_mean]                 , ...
                           'kwm', [   verify_kr(3),   kr_mean]);     

end
s = linspace(0, 1, 1001).'; kr = fluid.relperm(s);
plot(s, kr(:,1),'r'); hold on
plot(s, kr(:,2),'r'); hold on
