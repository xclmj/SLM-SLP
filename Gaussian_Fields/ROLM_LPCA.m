%% Gradient-based History Matching by finite difference
Jflux=0;
Jwatercut=0;
Jpres=0;

clear Jtotal_Outer_L_SeismicData_ProductionData Jtotal_L_OLSPCA_SeismicData_ProductionData update_KL_SeismicData_ProductionData
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

Vandermonde_Chebyshev=3;    % 1 means Vandermonde perturbation, 2 means Chebyshev perturbation, 3 means our perturbation
DP_IP=1;    %  1 means dependent perturbation, 2 means independent perturbation
FOM_ROM=2;  % 1 means FOM, 2 means ROM for cost function calculation
OPtimization_Algorithm=1;  % 1 means steepset descent, 2 means BFGS, 3 means Gaussian_Newton

Domain_Dx=[3 4 5];    Domain_Dy=[4 5 6];  

Num_kr=6;     % updating 6 relative permeability parameters
Both_K_Kr=2;  % 1 only updating log-permeability, 2 simultaneously updating 6 relative permeability parameters

for IIJ=1:1
    
    Jtotal_G=[];
    
    
    AAA=['------ Multi-lelvel:' num2str(IIJ)  '------'];
    disp(AAA)       

    Domain_Numx= Domain_Dx(IIJ);
    Domain_Numy= Domain_Dy(IIJ); 

    Initilization_DD1;
    
    Num_LPCA=(ceil(size(Layer{z}.fala,2)/Domain_Numx/Domain_Numy))*ones(Domain_Numx*Domain_Numy,1);
    
    parameterization_new1;
    Perturbation_Design_new;

   Num_K=0;
   for z=1:nz
        Num_K=Num_K+sum(Layer{z}.falaK);
   end

    verify_Linearization = zeros(Num_K+Num_kr,1); 
     verify_kr_initial=rand(Num_kr,1);
%     verify_kr_initial=ones(Num_kr,1);
    verify_Linearization(Num_K+1:Num_K+Num_kr,1)=verify_kr_initial;

    Perturb_Matrix_kr=delta*eye(Num_kr,Num_kr);
    NN=0;

    for ROLM=1:10

        AAA=['------ Outer-loop:' num2str(ROLM)  '------'];
        disp(AAA)  
        
        num_peturb=size(Perturb_Matrix,1)+size(Perturb_Matrix_kr,1)+1;
        % constructing linear model
        K_Layer=[];
        for jj=1:num_peturb

            if jj<size(Perturb_Matrix,1)+1
               verify2=verify_Linearization(1:Num_K)+Perturb_Matrix(jj,:)';
            else
               verify2=verify_Linearization(1:Num_K);
            end   

            for z=1:nz
                xx=0;
                for zz=1:z
                    xx=xx+Layer{zz}.D;
                end
                verify1=verify2(xx-Layer{zz}.D+1:xx);
                K=Layer{z}.mPermeablility+Layer{z}.fala*Layer{z}.TMGL*verify1;
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

        Gra_Flux=zeros(nWP,Num_K+Num_kr,ddd_welldata);
        for j=1:nWP   % for well-level
            for i=1:ddd_welldata
                B=zeros(Num_K,1);
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
                
                AA=Perturb_Matrix(:,BB);
                b=zeros(length(AA(1,:)),1);
                A=zeros(length(AA(1,:)),length(AA(1,:)));
                for k=1:size(Perturb_Matrix,1)
                    b=b+(FluxFD(j+(k-1)*nWP,i)-FluxFD(j+(num_peturb-1)*nWP,i))*AA(k,:)';
                    A=A+AA(k,:)'*AA(k,:);
                end 
                Gra_Flux(j,BB,i)=pinv(A)*b;
                
                for k=size(Perturb_Matrix,1)+1:num_peturb-1
                    Gra_Flux(j,k-size(Perturb_Matrix,1)+Num_K,i)= (FluxFD(j+(k-1)*nWP,i)-FluxFD(j+(num_peturb-1)*nWP,i))/0.01; 
                end
                
            end
        end
        
        Gra_Wat=zeros(nWP,Num_K+Num_kr,ddd_welldata);
        for j=1:nWP   % for well-level
            for i=1:ddd_welldata
                B=zeros(Num_K,1);
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
                AA=Perturb_Matrix(:,BB);
                b=zeros(length(AA(1,:)),1);
                A=zeros(length(AA(1,:)),length(AA(1,:)));
                for k=1:size(Perturb_Matrix,1)
                    b=b+(WatercutFD(j+(k-1)*nWP,i)-WatercutFD(j+(num_peturb-1)*nWP,i))*AA(k,:)';
                    A=A+AA(k,:)'*AA(k,:);
                end 
                Gra_Wat(j,BB,i)=pinv(A)*b;
                
                for k=size(Perturb_Matrix,1)+1:num_peturb-1
                    Gra_Wat(j,k-size(Perturb_Matrix,1)+Num_K,i)= (WatercutFD(j+(k-1)*nWP,i)-WatercutFD(j+(num_peturb-1)*nWP,i))/0.01; 
                end
                
            end
        end
        
        Gra_BHP=zeros(nWI,Num_K+Num_kr,ddd_welldata);
        for j=1:nWI   % for well-level
            for i=1:ddd_welldata
                B=zeros(Num_K,1);
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
                AA=Perturb_Matrix(:,BB);
                b=zeros(length(AA(1,:)),1);
                A=zeros(length(AA(1,:)),length(AA(1,:)));
                for k=1:size(Perturb_Matrix,1)
                    b=b+(BHPFD(j+(k-1)*nWI,i)-BHPFD(j+(num_peturb-1)*nWI,i))*AA(k,:)';
                    A=A+AA(k,:)'*AA(k,:);
                end 
                Gra_BHP(j,BB,i)=pinv(A)*b;
                
                for k=size(Perturb_Matrix,1)+1:num_peturb-1
                    Gra_BHP(j,k-size(Perturb_Matrix,1)+Num_K,i)= (BHPFD(j+(k-1)*nWI,i)-BHPFD(j+(num_peturb-1)*nWI,i))/0.01; 
                end  
                
            end
        end      
        
        Gra_Seismic=zeros(nx*ny*nz,Num_K+Num_kr,ddd_seismicdata);
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
                B=zeros(Num_K,1);
                % extract the local full column rank perturbation matrix
                [BB]=RBF(j,Layer{1}.falaK,Domain_index,Domain_Numx,Domain_Numy,nz);
                AA=Perturb_Matrix(:,BB);
                for z=1:length(Position)
                    b=zeros(length(AA(1,:)),1);
                    A=zeros(length(AA(1,:)),length(AA(1,:)));
                    for k=1:size(Perturb_Matrix,1)
                        b=b+(SatFD(Position(z)+(k-1)*nx*ny*nz,i)-SatFD(Position(z)+(num_peturb-1)*nx*ny*nz,i))*AA(k,:)';
                        A=A+AA(k,:)'*AA(k,:);
                    end
                    Gra_Seismic(Position(z),BB,i)=pinv(A)*b;
                end
            end
        end
        for j=1:nx*ny*nz   % for well-level
            for i=1:ddd_seismicdata
                for k=size(Perturb_Matrix,1)+1:num_peturb-1
                    Gra_Seismic(j,k-size(Perturb_Matrix,1)+Num_K,i)= (SatFD(j+(k-1)*nx*ny*nz,i)-SatFD(j+(num_peturb-1)*nx*ny*nz,i))/0.01; 
                end  
            end
        end
        
        
        num_peturb=1;
        if ROLM==1
           verify=verify_Linearization;
%            verify=(rand(Num_K,1)-0.5*ones(Num_K,1))*1.5;
           verify_init=verify; 
           
            xyz=FOM_ROM;
            FOM_ROM=1;

            for z=1:nz
                xx=0;
                for zz=1:z
                    xx=xx+Layer{zz}.D;
                end
                verify1=verify(xx-Layer{zz}.D+1:xx);
                K_upt=Layer{z}.mPermeablility+Layer{z}.fala*Layer{z}.TMGL*verify1;
                PORO=0.25 * (exp(K_upt) ./200).^0.1;
                ROCK=convertFrom(exp(K_upt),milli*darcy);
                rock{1}.perm(1+(z-1)*nx*ny:z*nx*ny)=ROCK;
                rock{1}.poro(1+(z-1)*nx*ny:z*nx*ny)=PORO;    

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
            
            Jtotal_Outer_L_SeismicData_ProductionData{IIJ}.CF(1,xyz)=Jtotal+Jseismic;
            Jtotal_Outer_L_SeismicData_ProductionData{IIJ}.Index(1,xyz)=1;
            Jtotal+Jseismic
            
            FOM_ROM=xyz;

        end

        Stop_Condition=0;
        for JJ=1:IteNum  %IteNum
            
            tic 
            
            verify2=verify;
            for z=1:nz
                xx=0;
                for zz=1:z
                    xx=xx+Layer{zz}.D;
                end
                verify1=verify2(xx-Layer{zz}.D+1:xx);
                K=Layer{z}.mPermeablility+Layer{z}.fala*Layer{z}.TMGL*verify1;
                PORO=0.25 * (exp(K) ./200).^0.1;
                ROCK=convertFrom(exp(K),milli*darcy);
                rock{1}.perm(1+(z-1)*nx*ny:z*nx*ny)=ROCK;
                rock{1}.poro(1+(z-1)*nx*ny:z*nx*ny)=PORO;    
            end

            if JJ==1 && ROLM==1
               K_init=log(convertTo(rock{1}.perm,milli*darcy));
            end
            
            Simulation_FOM_ROM_Seismic_ProductionData;
 
            % global cost function
            JtotalFDm_G=[];
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
                JtotalFDm_G=[JtotalFDm_G Jtotal];
                Jflux=0;
                Jwatercut=0;
                Jpres=0;
            end
            
            Jseismic=0;
            for tt=1:ddd_seismicdata
                tttt1=(1+(tt-1)*nx*ny*nz):(tt*nx*ny*nz);
                Jseismic=Jseismic+(SatFD(:,tt)-Sat(:,tt))'*inv(Rs(tttt1,tttt1))*(SatFD(:,tt)-Sat(:,tt))/2;
            end
            
            Jtotal_G=[Jtotal_G JtotalFDm_G(end)+Jseismic];
            JtotalFDm_G(end)+Jseismic
            
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

            Gradient_Computation_Seismic_ProductionData;
            
            if OPtimization_Algorithm<3
           
                if JJ==1 
                   ll=0.2;
                end

               if  JJ>1 && Jtotal_G(end)>Jtotal_G(end-1)
                    ll=ll/2;
                    verify=verifym-ll*FF;
                else 

                     verifym=verify;
                     verify=verifym-ll*EE;               
               end
                
            elseif OPtimization_Algorithm == 3
                
                   verifym=verify;
                   verify=verifym-EE;  
                
            end

            FF=EE;
            verify_upt=verify;
            AAA=['Inner-loop Iteration Step:' num2str(JJ)];
            disp(AAA)

            toc
            
        end
        Stop_Condition=0;

        for z=1:nz
            xx=0;
            for zz=1:z
                xx=xx+Layer{zz}.D;
            end
            verify1=verify(xx-Layer{zz}.D+1:xx);
            K_upt=Layer{z}.mPermeablility+Layer{z}.fala*Layer{z}.TMGL*verify1;
            PORO=0.25 * (exp(K_upt) ./200).^0.1;
            ROCK=convertFrom(exp(K_upt),milli*darcy);
            rock{1}.perm(1+(z-1)*nx*ny:z*nx*ny)=ROCK;
            rock{1}.poro(1+(z-1)*nx*ny:z*nx*ny)=PORO;    

        end        
        update_KL_SeismicData_ProductionData{IIJ}.Perm(1:length(rock{1}.perm),ROLM)= log(convertTo(rock{1}.perm,milli*darcy));
        update_KL_SeismicData_ProductionData{IIJ}.Re_Perm(1:Num_kr,ROLM)= verify(Num_K+1:Num_K+Num_kr,1);
        
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
               
        Jtotal_Outer_L_SeismicData_ProductionData{IIJ}.CF(ROLM+1,xyz)=Jtotal+Jseismic;
        Jtotal_Outer_L_SeismicData_ProductionData{IIJ}.Index(ROLM+1,xyz)=length(Jtotal_G);
        Jtotal+Jseismic


        FOM_ROM=xyz;


    end
    Jtotal_L_OLSPCA_SeismicData_ProductionData{IIJ}.G(1:length(Jtotal_G),FOM_ROM)=Jtotal_G';

    figure(1)
    subplot(1,2,1)
    semilogy(Jtotal_G); hold on
    if IIJ==1
        semilogy(1:IteNum*ROLM,ones(IteNum*ROLM,1)*((2*nWP+nWI)*ddd_welldata+size(SatFD,1)*size(SatFD,2))*5,'k--','LineWidth',2);
        semilogy(1:IteNum*ROLM,ones(IteNum*ROLM,1)*Jref_L,'r--','LineWidth',2);
    end
    
    subplot(1,2,2)
    semilogy(Jtotal_Outer_L_SeismicData_ProductionData{IIJ}.CF(:,xyz),'-o','LineWidth',1); hold on
    if IIJ==1
        semilogy(1:ROLM+1,ones(1+ROLM,1)*((2*nWP+nWI)*ddd_welldata+size(SatFD,1)*size(SatFD,2))*5,'k--','LineWidth',2); 
        semilogy(1:1+ROLM,ones(1+ROLM,1)*Jref_L,'r--','LineWidth',2);
    end
        
end


for i=1:nz
    
    K=K_ref(1+(i-1)*nx*ny:i*nx*ny,1);
    ROCK=nan*ones(nx*ny,1);
    ROCK(indexMap)= K(indexMap);
    ROCK=reshape(ROCK,nx,ny,1);
    figure('color',[1,1,1])
    cm = colormap(jet);
    gca=pcolor(ROCK); colorbar; box on;caxis([1, 8])
    set(gca, 'LineStyle','none');
    hold on
    plot(jInj, iInj, '^k', 'MarkerFaceColor', [0.99 0.99 0.99], 'MarkerSize', 8); % used to differentiate injection and production wells
    for index=1:length(jInj)
        text(jInj(index)-1, iInj(index)-2, ['I' num2str(index)], 'FontName', 'Times New Roman', 'color', 'k', 'FontSize', 12, 'FontWeight', 'b');
    end

    plot(jProd, iProd, 'ok', 'MarkerFaceColor',[0.99 0.99 0.99], 'MarkerSize',8);
    for index=1:length(jProd)
        text(jProd(index)-2, iProd(index)-2, ['P' num2str(index)], 'FontName', 'Times New Roman', 'color', 'k', 'FontSize', 12, 'FontWeight', 'b'); 
    end
    title([num2str(i) 'st layer True'])

    K=update_KL_SeismicData_ProductionData{IIJ}.Perm(1+(i-1)*nx*ny:i*nx*ny,ROLM);
    ROCK=nan*ones(nx*ny,1);
    ROCK(indexMap)= K(indexMap);
    ROCK=reshape(ROCK,nx,ny,1);
    figure('color',[1,1,1])
    cm = colormap(jet);
    gca=pcolor(ROCK); colorbar; box on;caxis([1, 8])
    set(gca, 'LineStyle','none');
    hold on
    plot(jInj, iInj, '^k', 'MarkerFaceColor', [0.99 0.99 0.99], 'MarkerSize', 8); % used to differentiate injection and production wells
    for index=1:length(jInj)
        text(jInj(index)-1, iInj(index)-2, ['I' num2str(index)], 'FontName', 'Times New Roman', 'color', 'k', 'FontSize', 12, 'FontWeight', 'b');
    end

    plot(jProd, iProd, 'ok', 'MarkerFaceColor',[0.99 0.99 0.99], 'MarkerSize',8);
    for index=1:length(jProd)
        text(jProd(index)-2, iProd(index)-2, ['P' num2str(index)], 'FontName', 'Times New Roman', 'color', 'k', 'FontSize', 12, 'FontWeight', 'b'); 
    end
    title([num2str(i) 'st layer Updated'])

    K=K_init(1+(i-1)*nx*ny:i*nx*ny,1);
    ROCK=nan*ones(nx*ny,1);
    ROCK(indexMap)= K(indexMap);
    ROCK=reshape(ROCK,nx,ny,1);
    figure('color',[1,1,1])
    cm = colormap(jet);
    gca=pcolor(ROCK); colorbar; box on;caxis([1, 8])
    set(gca, 'LineStyle','none');
    hold on
    plot(jInj, iInj, '^k', 'MarkerFaceColor', [0.99 0.99 0.99], 'MarkerSize', 8); % used to differentiate injection and production wells
    for index=1:length(jInj)
        text(jInj(index)-1, iInj(index)-2, ['I' num2str(index)], 'FontName', 'Times New Roman', 'color', 'k', 'FontSize', 12, 'FontWeight', 'b');
    end

    plot(jProd, iProd, 'ok', 'MarkerFaceColor',[0.99 0.99 0.99], 'MarkerSize',8);
    for index=1:length(jProd)
        text(jProd(index)-2, iProd(index)-2, ['P' num2str(index)], 'FontName', 'Times New Roman', 'color', 'k', 'FontSize', 12, 'FontWeight', 'b'); 
    end
    title([num2str(i) 'st layer Initial'])

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

%     verify_kr=(exp(verify_kr).*[n_max; n_max; S_max; S_max; kr_max; kr_max]+[n_min; n_min; S_min; S_min; kr_min; kr_min])./(ones(Num_kr,1)+exp(verify_kr));    
    verify_kr=verify_kr.*([n_max; n_max; S_max; S_max; kr_max; kr_max]-[n_min; n_min; S_min; S_min; kr_min; kr_min])+[n_min; n_min; S_min; S_min; kr_min; kr_min];
    fluid = initCoreyFluid('mu' , [ 0.4, 2]*centi*poise     , ...
                           'rho', [1014, 859]*kilogram/meter^3, ...
                           'n'  , [   verify_kr(1),   verify_kr(2)]                 , ...
                           'sr' , [ verify_kr(3), verify_kr(4)]                 , ...
                           'kwm', [  verify_kr(5),   verify_kr(6)]);

elseif Num_kr==3 && Both_K_Kr==2

    verify_kr=verify_kr_initial; 
%     verify_kr=(exp(verify_kr).*[n_max; n_max; kr_max]+[n_min; n_min; kr_min])./(ones(Num_kr,1)+exp(verify_kr));
    verify_kr=verify_kr.*([n_max; n_max; S_max; S_max; kr_max; kr_max]-[n_min; n_min; S_min; S_min; kr_min; kr_min])+[n_min; n_min; S_min; S_min; kr_min; kr_min];
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

%     verify_kr=(exp(verify_kr).*[n_max; n_max; S_max; S_max; kr_max; kr_max]+[n_min; n_min; S_min; S_min; kr_min; kr_min])./(ones(Num_kr,1)+exp(verify_kr));    
    verify_kr=verify_kr.*([n_max; n_max; S_max; S_max; kr_max; kr_max]-[n_min; n_min; S_min; S_min; kr_min; kr_min])+[n_min; n_min; S_min; S_min; kr_min; kr_min];

    fluid = initCoreyFluid('mu' , [ 0.4, 2]*centi*poise     , ...
                           'rho', [1014, 859]*kilogram/meter^3, ...
                           'n'  , [   verify_kr(1),   verify_kr(2)]                 , ...
                           'sr' , [ verify_kr(3), verify_kr(4)]                 , ...
                           'kwm', [  verify_kr(5),   verify_kr(6)]);

elseif Num_kr==3 && Both_K_Kr==2

    verify_kr=verify(Num_K+1:Num_K+Num_kr,1); 
%     verify_kr=(exp(verify_kr).*[n_max; n_max; kr_max]+[n_min; n_min; kr_min])./(ones(Num_kr,1)+exp(verify_kr));
    verify_kr=verify_kr.*([n_max; n_max; S_max; S_max; kr_max; kr_max]-[n_min; n_min; S_min; S_min; kr_min; kr_min])+[n_min; n_min; S_min; S_min; kr_min; kr_min];

    fluid = initCoreyFluid('mu' , [ 0.4, 2]*centi*poise     , ...
                           'rho', [1014, 859]*kilogram/meter^3, ...
                           'n'  , [   verify_kr(1),   verify_kr(2)]                 , ...
                           'sr' , [ S_mean, S_mean]                 , ...
                           'kwm', [   verify_kr(3),   kr_mean]);     

end
s = linspace(0, 1, 1001).'; kr = fluid.relperm(s);
plot(s, kr(:,1),'r'); hold on
plot(s, kr(:,2),'r'); hold on
