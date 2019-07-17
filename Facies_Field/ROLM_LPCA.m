
disp('SLM-SLP without adaptive scheme')
%% Gradient-based History Matching by finite difference
Jflux=0;
Jwatercut=0;
Jpres=0;

clear Jtotal_Outer_L Jtotal_L_OLSPCA update_KL
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

Vandermonde_Chebyshev=3;    % 1 means Vandermonde perturbation, 2 means Chebyshev perturbation, 3 means our perturbation
DP_IP=1;    %  1 means dependent perturbation, 2 means independent perturbation
FOM_ROM=2;  % 1 means FOM, 2 means ROM for cost function calculation
OPtimization_Algorithm=1;  % 1 means steepset descent, 2 means BFGS, 3 means Gaussian_Newton
Domain_Dx=[3 4 5 2 8];    Domain_Dy=[3 4 5 8 2];  

for IIJ=1:1
    
    Jtotal_G=[];
       
    AAA=['------ Domain Decomposition:' num2str(Domain_Dx(IIJ)) '\times' num2str(Domain_Dy(IIJ)) '------'];
    disp(AAA)       

    Domain_Numx=Domain_Dx(IIJ);
    Domain_Numy=Domain_Dy(IIJ); 

    Initilization_DD1;
    
    EnK_num=(ceil(size(Layer{z}.fala,2)/Domain_Numx/Domain_Numy));

    parameterization_new1;
    Perturbation_Design_new;

    Num_K=0;
    for z=1:nz
        Num_K=Num_K+sum(Layer{z}.falaK);
    end

    verify_Linearization=zeros(Num_K,1);  

    num_peturb=size(Perturb_Matrix,1)+1;
    NN=0;

    for ROLM=1:10

        AAA=['------ Outer-loop:' num2str(ROLM)  '------'];
        disp(AAA)  
        
        num_peturb=size(Perturb_Matrix,1)+1;
        % constructing linear model
        K_Layer=[];
        for jj=1:num_peturb

            if jj<num_peturb
               verify2=verify_Linearization+Perturb_Matrix(jj,:)';
            else
               verify2=verify_Linearization;
            end   

            for z=1:nz
                xx=0;
                for zz=1:z
                    xx=xx+Layer{zz}.D;
                end
                verify1=verify2(xx-Layer{zz}.D+1:xx);
                K=Layer{z}.mPermeablility+Layer{z}.fala*Layer{z}.TMGL*verify1;
                
                [K,S]= OPCA_TwoFacies(Layer{z}.fala, Layer{z}.TMGL*verify1, Layer{z}.mPermeablility, gammaF);              

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

        Gra_Flux=zeros(nWP,Num_K,ddd_welldata);
        for j=1:nWP   % for well-level
            for i=1:ddd_welldata
                B=zeros(Num_K,1);
                % extract the local full column rank perturbation matrix
                [BB]=RBF(WellP_Domain(j),Layer{1}.falaK,Domain_index,Domain_Numx,Domain_Numy,nz);
                AA=Perturb_Matrix(:,BB);
                b=zeros(length(AA(1,:)),1);
                A=zeros(length(AA(1,:)),length(AA(1,:)));
                for k=1:num_peturb-1
                    b=b+(FluxFD(j+(k-1)*nWP,i)-FluxFD(j+(num_peturb-1)*nWP,i))*AA(k,:)';
                    A=A+AA(k,:)'*AA(k,:);
                end 
                Gra_Flux(j,BB,i)=pinv(A)*b;
            end
        end
        
        Gra_Wat=zeros(nWP,Num_K,ddd_welldata);
        for j=1:nWP   % for well-level
            for i=1:ddd_welldata
                B=zeros(Num_K,1);
                % extract the local full column rank perturbation matrix
                [BB]=RBF(WellP_Domain(j),Layer{1}.falaK,Domain_index,Domain_Numx,Domain_Numy,nz);
                AA=Perturb_Matrix(:,BB);
                b=zeros(length(AA(1,:)),1);
                A=zeros(length(AA(1,:)),length(AA(1,:)));
                for k=1:num_peturb-1
                    b=b+(WatercutFD(j+(k-1)*nWP,i)-WatercutFD(j+(num_peturb-1)*nWP,i))*AA(k,:)';
                    A=A+AA(k,:)'*AA(k,:);
                end 
                Gra_Wat(j,BB,i)=pinv(A)*b;
            end
        end
        
        Gra_BHP=zeros(nWI,Num_K,ddd_welldata);
        for j=1:nWI   % for well-level
            for i=1:ddd_welldata
                B=zeros(Num_K,1);
                % extract the local full column rank perturbation matrix
                [BB]=RBF(WellI_Domain(j),Layer{1}.falaK,Domain_index,Domain_Numx,Domain_Numy,nz);
                AA=Perturb_Matrix(:,BB);
                b=zeros(length(AA(1,:)),1);
                A=zeros(length(AA(1,:)),length(AA(1,:)));
                for k=1:num_peturb-1
                    b=b+(BHPFD(j+(k-1)*nWI,i)-BHPFD(j+(num_peturb-1)*nWI,i))*AA(k,:)';
                    A=A+AA(k,:)'*AA(k,:);
                end 
                Gra_BHP(j,BB,i)=pinv(A)*b;
            end
        end        
        
        
        num_peturb=1;
        if ROLM==1
           verify=zeros(Num_K,1);
           verify_init=verify; 
           
            xyz=FOM_ROM;
            FOM_ROM=1;

            for z=1:nz
                xx=0;
                for zz=1:z
                    xx=xx+Layer{zz}.D;
                end
                verify1=verify(xx-Layer{zz}.D+1:xx);
                K=Layer{z}.mPermeablility+Layer{z}.fala*Layer{z}.TMGL*verify1;
                
                if Facies_tyep<4
                   [K,S]= OPCA_TwoFacies(Layer{z}.fala, Layer{z}.TMGL*verify1, Layer{z}.mPermeablility, gammaF);
                else
                   [K,S] = OPCA_ThreeFacies(Layer{z}.fala, Layer{z}.TMGL*verify1, Layer{z}.mPermeablility, gamma11, gamma12, gamma21, gamma22);
                end 

                K=K_mud*exp(log(K_sand/K_mud)*K); 
                PORO=0.25 * (K ./200).^0.1;
                ROCK=convertFrom(K,milli*darcy);
                rock{1}.perm(1+(z-1)*nx*ny:z*nx*ny)=ROCK;
                rock{1}.poro(1+(z-1)*nx*ny:z*nx*ny)=PORO;    

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
            Jtotal_Outer_L{IIJ}.CF(1,xyz)=Jtotal;
            Jtotal_Outer_L{IIJ}.Index(1,xyz)=1;
            Jtotal
            
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
                [K,S]= OPCA_TwoFacies(Layer{z}.fala, Layer{z}.TMGL*verify1, Layer{z}.mPermeablility, gammaF);
                K1=K_mud*exp(log(K_sand/K_mud)*K); 
                PORO=0.25 * (K1 ./200).^0.1;
                ROCK=convertFrom(K1,milli*darcy);
                rock{1}.perm(1+(z-1)*nx*ny:z*nx*ny)=ROCK;
                rock{1}.poro(1+(z-1)*nx*ny:z*nx*ny)=PORO;    
            end

            if JJ==1 && ROLM==1
               K_init=K;
            end
            Simulation_FOM_ROM;
 
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
            Jtotal_G=[Jtotal_G JtotalFDm_G(end)];
            JtotalFDm_G(end)
            
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

            
            Gradient_Computation;
            
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
            [K_upt,S]= OPCA_TwoFacies(Layer{z}.fala, Layer{z}.TMGL*verify1, Layer{z}.mPermeablility, gammaF);
            K=K_mud*exp(log(K_sand/K_mud)*K_upt); 
            PORO=0.25 * (K ./200).^0.1;
            ROCK=convertFrom(K,milli*darcy);
            rock{1}.perm(1+(z-1)*nx*ny:z*nx*ny)=ROCK;
            rock{1}.poro(1+(z-1)*nx*ny:z*nx*ny)=PORO;    

        end        
        update_KL{IIJ}.Perm(1:length(rock{1}.perm),ROLM)= K_upt;
        
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
        Jtotal_Outer_L{IIJ}.CF(ROLM+1,xyz)=Jtotal;
        Jtotal_Outer_L{IIJ}.Index(ROLM+1,xyz)=length(Jtotal_G);
        Jtotal

        FOM_ROM=xyz;

    end
    Jtotal_L_OLSPCA{IIJ}.G(1:length(Jtotal_G),FOM_ROM)=Jtotal_G';

    figure(1)
    subplot(1,2,1)
    semilogy(Jtotal_G); hold on
    if IIJ==1
        semilogy(1:IteNum*ROLM,ones(IteNum*ROLM,1)*(2*nWP+nWI)*ddd_welldata*5/2,'k--','LineWidth',2);
        semilogy(1:IteNum*ROLM,ones(IteNum*ROLM,1)*Jref_L,'r--','LineWidth',2);
    end
    
    subplot(1,2,2)
    semilogy(Jtotal_Outer_L{IIJ}.CF(:,xyz),'-o','LineWidth',1); hold on
    if IIJ==1
        semilogy(1:ROLM+1,ones(1+ROLM,1)*(2*nWP+nWI)*ddd_welldata*5/2,'k--','LineWidth',2); 
        semilogy(1:1+ROLM,ones(1+ROLM,1)*Jref_L,'r--','LineWidth',2);
    end
        
end


for i=1:nz
    
    K=K_ref(1+(i-1)*nx*ny:i*nx*ny,1);
    ROCK=nan*ones(nx*ny,1);
    ROCK= K;
    ROCK=reshape(ROCK,nx,ny,1);
    figure('color',[1,1,1])
    cm = colormap(jet);
    gca=pcolor(ROCK); colorbar; box on;caxis([0, 1])
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

    K=update_KL{IIJ}.Perm(1+(i-1)*nx*ny:i*nx*ny,ROLM);
    ROCK=nan*ones(nx*ny,1);
    ROCK= K;
    ROCK=reshape(ROCK,nx,ny,1);
    figure('color',[1,1,1])
    cm = colormap(jet);
    gca=pcolor(ROCK); colorbar; box on;caxis([0, 1])
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
    ROCK= K;
    ROCK=reshape(ROCK,nx,ny,1);
    figure('color',[1,1,1])
    cm = colormap(jet);
    gca=pcolor(ROCK); colorbar; box on;caxis([0, 1])
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
