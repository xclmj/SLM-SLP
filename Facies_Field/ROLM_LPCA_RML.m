%% Gradient-based History Matching by finite difference
Jflux=0;
Jwatercut=0;
Jpres=0;
  
clear Jtotal_Outer_RML Jtotal_RML_OLSPCA update_KRML
% global cost function at the reference field
for ii=1:Num_RML
    Jflux=0;
    Jwatercut=0;
    Jpres=0;
    tttt3=(1+(ii-1)*nWI):(ii*nWI);
    tttt4=(1+(ii-1)*nWP):(ii*nWP);
    for tt=1:ddd_welldata
        tttt1=(1+(tt-1)*nWI):(nWI+(tt-1)*nWI);
        Jpres=Jpres+(BHPref(:,tt)-BHP_RML(tttt3,tt))'*inv(Rp(tttt1,tttt1))*(BHPref(:,tt)-BHP_RML(tttt3,tt))/2;
        tttt2=(1+(tt-1)*nWP):(nWP+(tt-1)*nWP);  
        Jwatercut=Jwatercut+(Watercutref(:,tt)-Watercut_RML(tttt4,tt))'*inv(Rfw(tttt2,tttt2))*(Watercutref(:,tt)-Watercut_RML(tttt4,tt))/2; 
        Jflux=Jflux+(Fluxref(:,tt)-Flux_RML(tttt4,tt))'*inv(Rq(tttt2,tttt2))*(Fluxref(:,tt)-Flux_RML(tttt4,tt))/2;
    end    
    Jtotal=Jflux+Jwatercut+Jpres;
    Jref_L_RML(ii)= Jtotal;
end

% local cost function at the reference field
erro1=1e-4;  % iteration criteration
erro2=1e-3;
IteNum=100;
sigma1=1;
sigma2=1;

Vandermonde_Chebyshev=3;    % 1 means Vandermonde perturbation, 2 means Chebyshev perturbation
DP_IP=1;    %  1 means dependent perturbation, 2 means independent perturbation
FOM_ROM=2;  % 1 means FOM, 2 means ROM for cost function calculation
OPtimization_Algorithm=1;  % 1 means steepset descent, 2 means BFGS, 3 means Gaussian_Newton
XXX=[2 5 8];

Domain_Numx=3;
Domain_Numy=3;
EnK_num=5;

Initilization_DD1;
parameterization_new1;
Perturbation_Design_new;
Num_K=0;
for z=1:nz
    Num_K=Num_K+sum(Layer{z}.falaK);
end

for ii=1:Num_RML
    
    verify_RML_Init(1:Num_K,ii)=(rand(Num_K,1)-0.5*ones(Num_K,1))*2;
    
end

verify_RML=verify_RML_Init;
for IIJ=1:Num_RML
    
    Jtotal_G=[];
    
    
    AAA=['------ RML:' num2str(IIJ)  '------'];
    disp(AAA)       

%     verify_Linearization=zeros(Num_K,1);  
    verify_Linearization=verify_RML(1:Num_K,IIJ);

    num_peturb=size(Perturb_Matrix,1)+1;

    for ROLM=1:10

        AAA=['------ Outer-loop:' num2str(ROLM)  '------'];
        disp(AAA)  
        
        num_peturb=size(Perturb_Matrix,1)+1;
        % constructing linear model
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

            verify=verify_RML(:,IIJ);
           
            xyz=FOM_ROM;
            FOM_ROM=1;

            for z=1:nz
                xx=0;
                for zz=1:z
                    xx=xx+Layer{zz}.D;
                end
                verify1=verify(xx-Layer{zz}.D+1:xx);
                K=Layer{z}.mPermeablility+Layer{z}.fala*Layer{z}.TMGL*verify1;

                [K,S]= OPCA_TwoFacies(Layer{z}.fala, Layer{z}.TMGL*verify1, Layer{z}.mPermeablility, gammaF);
                
                K=K_mud*exp(log(K_sand/K_mud)*K); 
                PORO=0.25 * (K ./200).^0.1;
                ROCK=convertFrom(K,milli*darcy);
                rock{1}.perm(1+(z-1)*nx*ny:z*nx*ny)=ROCK;
                rock{1}.poro(1+(z-1)*nx*ny:z*nx*ny)=PORO;   

            end

            Simulation_FOM_ROM;

            Jflux=0;
            Jwatercut=0;
            Jpres=0;
            tttt1=(1+(IIJ-1)*nWI):(nWI+(IIJ-1)*nWI);
            tttt2=(1+(IIJ-1)*nWP):(length(rSol.wellSol)-nWI+(IIJ-1)*nWP);
            for tt=1:ddd_welldata
                tttt3=(1+(tt-1)*nWI):(nWI+(tt-1)*nWI);
                Jpres=Jpres+(BHP_RML(tttt1,tt)-BHPFD(:,tt))'*inv(Rp(tttt3,tttt3))*(BHP_RML(tttt1,tt)-BHPFD(:,tt))/2;
                tttt4=(1+(tt-1)*nWP):(nWP+(tt-1)*nWP);
                Jflux=Jflux+(Flux_RML(tttt2,tt)-FluxFD(:,tt))'*inv(Rq(tttt4,tttt4))*(Flux_RML(tttt2,tt)-FluxFD(:,tt))/2;
                Jwatercut=Jwatercut+(Watercut_RML(tttt2,tt)-WatercutFD(:,tt))'*inv(Rfw(tttt4,tttt4))*(Watercut_RML(tttt2,tt)-WatercutFD(:,tt))/2;         
            end  
            Jtotal=Jflux+Jwatercut+Jpres;
            Jflux=0;
            Jwatercut=0;
            Jpres=0;

            Jtotal_Outer_RML{IIJ}.CF(1,xyz)=Jtotal;
            Jtotal_Outer_RML{IIJ}.Index(1,xyz)=1;
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
                
                K=K_mud*exp(log(K_sand/K_mud)*K); 
                PORO=0.25 * (K ./200).^0.1;
                ROCK=convertFrom(K,milli*darcy);
                rock{1}.perm(1+(z-1)*nx*ny:z*nx*ny)=ROCK;
                rock{1}.poro(1+(z-1)*nx*ny:z*nx*ny)=PORO;    
            end

            if JJ==1 && ROLM==1
               K_init_RML(:,IIJ)=log(convertTo(rock{1}.perm,milli*darcy));
            end
            Simulation_FOM_ROM;
 
            % global cost function
            JtotalFDm_G=[];

            Jflux=0;
            Jwatercut=0;
            Jpres=0;
            tttt1=(1+(IIJ-1)*nWI):(nWI+(IIJ-1)*nWI);
            tttt2=(1+(IIJ-1)*nWP):(length(rSol.wellSol)-nWI+(IIJ-1)*nWP);
            for tt=1:ddd_welldata
                tttt3=(1+(tt-1)*nWI):(nWI+(tt-1)*nWI);
                Jpres=Jpres+(BHP_RML(tttt1,tt)-BHPFD(:,tt))'*inv(Rp(tttt3,tttt3))*(BHP_RML(tttt1,tt)-BHPFD(:,tt))/2;
                tttt4=(1+(tt-1)*nWP):(nWP+(tt-1)*nWP);
                Jflux=Jflux+(Flux_RML(tttt2,tt)-FluxFD(:,tt))'*inv(Rq(tttt4,tttt4))*(Flux_RML(tttt2,tt)-FluxFD(:,tt))/2;
                Jwatercut=Jwatercut+(Watercut_RML(tttt2,tt)-WatercutFD(:,tt))'*inv(Rfw(tttt4,tttt4))*(Watercut_RML(tttt2,tt)-WatercutFD(:,tt))/2;         
            end  
            Jtotal=Jflux+Jwatercut+Jpres;
            Jflux=0;
            Jwatercut=0;
            Jpres=0;
            JtotalFDm_G=[JtotalFDm_G Jtotal];

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

            
            Gradient_Computation_RML;
            
            if OPtimization_Algorithm<3
           
                if JJ==1 
                   ll=0.1;
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
        update_KRML{IIJ}.Perm(1:length(K_upt),ROLM)= K_upt;
        
        verify_Linearization=verify;
        
        
        
        xyz=FOM_ROM;
        FOM_ROM=1;

        Simulation_FOM_ROM;

        Jflux=0;
        Jwatercut=0;
        Jpres=0;
        tttt1=(1+(IIJ-1)*nWI):(nWI+(IIJ-1)*nWI);
        tttt2=(1+(IIJ-1)*nWP):(length(rSol.wellSol)-nWI+(IIJ-1)*nWP);
        for tt=1:ddd_welldata
            tttt3=(1+(tt-1)*nWI):(nWI+(tt-1)*nWI);
            Jpres=Jpres+(BHP_RML(tttt1,tt)-BHPFD(:,tt))'*inv(Rp(tttt3,tttt3))*(BHP_RML(tttt1,tt)-BHPFD(:,tt))/2;
            tttt4=(1+(tt-1)*nWP):(nWP+(tt-1)*nWP);
            Jflux=Jflux+(Flux_RML(tttt2,tt)-FluxFD(:,tt))'*inv(Rq(tttt4,tttt4))*(Flux_RML(tttt2,tt)-FluxFD(:,tt))/2;
            Jwatercut=Jwatercut+(Watercut_RML(tttt2,tt)-WatercutFD(:,tt))'*inv(Rfw(tttt4,tttt4))*(Watercut_RML(tttt2,tt)-WatercutFD(:,tt))/2;         
        end  
        Jtotal=Jflux+Jwatercut+Jpres;
        Jflux=0;
        Jwatercut=0;
        Jpres=0;
        Jtotal_Outer_RML{IIJ}.CF(ROLM+1,xyz)=Jtotal;
        Jtotal_Outer_RML{IIJ}.Index(ROLM+1,xyz)=length(Jtotal_G);
        Jtotal


        FOM_ROM=xyz;


    end
     Jtotal_RML_OLSPCA{IIJ}.G(1:length(Jtotal_G),FOM_ROM)=Jtotal_G';
    
    figure(1)
    semilogy(Jtotal_Outer_RML{IIJ}.CF(:,xyz),'-o','LineWidth',1); hold on
        
end


