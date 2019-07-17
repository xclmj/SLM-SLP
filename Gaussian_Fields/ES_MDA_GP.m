%% Data Assimilation methods: EnRML, and ES-MDA
Domain_Numx=3;
Domain_Numy=3;
Initilization_DD1;
EnK_num=5;
Num_LPCA=EnK_num*ones(Domain_Numx*Domain_Numy,1);
parameterization_new1;

Data_type=2;   % 1 means well data, 2 means seismic data, 3 means both of well and seismic data

Num_RML=500;
%% Step 1: generating initial ensemble
% initial ensemble parameters
Num_K=0;
for z=1:nz
    Num_K=Num_K+sum(Layer{z}.falaK);
end

verify_RML_Init=[];
for ii=1:Num_RML
    
    verify_RML_Init(1:Num_K,ii)=(rand(Num_K,1)-0.5*ones(Num_K,1))*2;
    verify_RML_Init(Num_K+1:Num_K+Num_kr,ii)=rand(Num_kr,1);
    
end
verify_RML=verify_RML_Init;
   
%% ES-MDA
Iteration_ESMDA=10;
RMSE_K_ESMDA=[];
RMSE_ensemble_ESMDA=[];
RMSE_Kr1_ESMDA=[];
RMSE_Kr2_ESMDA=[];
RMSE_Kr3_ESMDA=[];
RMSE_Kr4_ESMDA=[];
RMSE_Kr5_ESMDA=[];
RMSE_Kr6_ESMDA=[];
ESP_K_ESMDA=[];
        

Num_RML=100;
verify_ESMDA=verify_RML_Init(:,1:Num_RML);
clear rock
for Ite=1:Iteration_ESMDA

   AAA=['ES-MDA Iteration Step:' num2str(Ite)];
   disp(AAA)

   %% STEP 1: running ensemble     
   for jj=1:Num_RML

        verify2=verify_ESMDA(1:Num_K,jj);
        verify_Linearization=verify_ESMDA(:,jj);

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

    num_peturb = Num_RML;
    
    tic
    Gradient_Calculation_finitedifferenceESMDA;
    toc

    Flux_Simulation=[];
    Watercut_Simulation=[];
    BHP_Simulation=[];
    Sat_Simulation=[];
    
    % ensemble of simulated measurements
    for jj=1:Num_RML
        
        % ensemble of simulated measurements
        Flux11 = FluxFD(1+(jj-1)*nWP:jj*nWP,:);
        Watercut11 = WatercutFD(1+(jj-1)*nWP:jj*nWP,:);
        Pres11 = BHPFD(1+(jj-1)*nWI:jj*nWI,:);
        Sat11 = SatFD(1+(jj-1)*nx*ny*nz:jj*nx*ny*nz,:);
        
        Flux_Simulation=[Flux_Simulation  reshape(Flux11,size(Flux11,1)*size(Flux11,2),1)];
        Watercut_Simulation=[Watercut_Simulation  reshape(Watercut11,size(Watercut11,1)*size(Watercut11,2),1)];
        BHP_Simulation=[BHP_Simulation  reshape(Pres11,size(Pres11,1)*size(Pres11,2),1)];  
        Sat_Simulation=[Sat_Simulation  reshape(Sat11,size(Sat11,1)*size(Sat11,2),1)];  

    end
    clear Flux11 Watercut11 Pres11
    

    % ensemble of perturbed measurement
    Flux_Perturbation1 = reshape(Fluxref,size(Fluxref,1)*size(Fluxref,2),1);
    Watercut_Perturbation1 = reshape(Watercutref,size(Watercutref,1)*size(Watercutref,2),1);
    BHP_Perturbation1 = reshape(BHPref,size(BHPref,1)*size(BHPref,2),1);  
    Sat_Perturbation1 = reshape(Satref,size(Satref,1)*size(Satref,2),1);  

    Flux_Perturbation=[];
    Watercut_Perturbation=[];
    BHP_Perturbation=[];
    Sat_Perturbation=[];
    for jj=1:Num_RML
        
        Flux_Perturbation=[Flux_Perturbation  normrnd(Flux_Perturbation1,Flux_Perturbation1*0.05*sqrt(Iteration_ESMDA),size(Fluxref,1)*size(Fluxref,2),1)];
        Watercut_Perturbation=[Watercut_Perturbation  normrnd(Watercut_Perturbation1, ones(size(Watercutref,1)*size(Watercutref,2),1)*0.05*sqrt(Iteration_ESMDA),size(Watercutref,1)*size(Watercutref,2),1)];
        BHP_Perturbation=[BHP_Perturbation normrnd(BHP_Perturbation1,BHP_Perturbation1*0.05*sqrt(Iteration_ESMDA),size(BHPref,1)*size(BHPref,2),1)];
        Sat_Perturbation=[Sat_Perturbation normrnd(Sat_Perturbation1,Sat_Perturbation1*0.05*sqrt(Iteration_ESMDA),size(Satref,1)*size(Satref,2),1)];
    end
    clear Flux_Perturbation1 Watercut_Perturbation1 BHP_Perturbation1 Sat_Perturbation1
    
    
    xx=0;
    for jj=1:Num_RML
        
        xx=0;
        if Data_type==1
            xx=xx+(Flux_Simulation(:,jj)-Flux_Perturbation(:,jj))'*inv(Rq)*(Flux_Simulation(:,jj)-Flux_Perturbation(:,jj));
            xx=xx+(Watercut_Simulation(:,jj)-Watercut_Perturbation(:,jj))'*inv(Rfw)*(Watercut_Simulation(:,jj)-Watercut_Perturbation(:,jj));
            xx=xx+(BHP_Simulation(:,jj)-BHP_Perturbation(:,jj))'*inv(Rp)*(BHP_Simulation(:,jj)-BHP_Perturbation(:,jj));
        elseif Data_type==2
            xx=(Sat_Simulation(:,jj)-Sat_Perturbation(:,jj))'*inv(Rs)*(Sat_Simulation(:,jj)-Sat_Perturbation(:,jj));
        else
            xx=xx+(Flux_Simulation(:,jj)-Flux_Perturbation(:,jj))'*inv(Rq)*(Flux_Simulation(:,jj)-Flux_Perturbation(:,jj));
            xx=xx+(Watercut_Simulation(:,jj)-Watercut_Perturbation(:,jj))'*inv(Rfw)*(Watercut_Simulation(:,jj)-Watercut_Perturbation(:,jj));
            xx=xx+(BHP_Simulation(:,jj)-BHP_Perturbation(:,jj))'*inv(Rp)*(BHP_Simulation(:,jj)-BHP_Perturbation(:,jj));
            xx=xx+(Sat_Simulation(:,jj)-Sat_Perturbation(:,jj))'*inv(Rs)*(Sat_Simulation(:,jj)-Sat_Perturbation(:,jj));
        end
        RMSE_ensemble_ESMDA(jj,Ite)=xx;
    end
    
    for jj=1:Num_RML
        RMSE_K_ESMDA(jj,Ite)=norm(log(convertTo(rock{jj}.perm,milli*darcy))-K_ref,2)/sqrt(length(K_ref));
    end 
    
   
    for jj=1:Num_RML
        
        verify_kr=verify_ESMDA(Num_K+1:Num_K+Num_kr,jj);
        verify_kr=verify_kr.*([n_max; n_max; S_max; S_max; kr_max; kr_max]-[n_min; n_min; S_min; S_min; kr_min; kr_min])+[n_min; n_min; S_min; S_min; kr_min; kr_min];
        RMSE_Kr1_ESMDA(jj,Ite)=verify_kr(1);
        RMSE_Kr2_ESMDA(jj,Ite)=verify_kr(2);
        RMSE_Kr3_ESMDA(jj,Ite)=verify_kr(3);
        RMSE_Kr4_ESMDA(jj,Ite)=verify_kr(4);
        RMSE_Kr5_ESMDA(jj,Ite)=verify_kr(5);
        RMSE_Kr6_ESMDA(jj,Ite)=verify_kr(6);
    end

    mean_k=[];
    for jj=1:Num_RML
        mean_k=[mean_k log(convertTo(rock{jj}.perm,milli*darcy))];
    end
    mean_kk=mean(mean_k,2);
    mean_k=mean_k-repmat(mean_kk,1,Num_RML);
    
    kkk=0;
    for ii=1:length(K_ref)
        kkk=kkk+norm(mean_k(ii,:),2)/sqrt(Num_RML);
    end
    ESP_K_ESMDA(Ite)=kkk/length(K_ref);
    

    
    
    Perturb_verify_ESMDA = (verify_ESMDA-repmat(mean(verify_ESMDA,2),1,Num_RML))/sqrt(Num_RML-1);

    Perturb_Flux_Simulation = (Flux_Simulation-repmat(mean(Flux_Simulation,2),1,Num_RML))/sqrt(Num_RML-1);
    Perturb_Watercut_Simulation = (Watercut_Simulation-repmat(mean(Watercut_Simulation,2),1,Num_RML))/sqrt(Num_RML-1);
    Perturb_BHP_Simulation = (BHP_Simulation-repmat(mean(BHP_Simulation,2),1,Num_RML))/sqrt(Num_RML-1);
    Perturb_Sat_Simulation = (Sat_Simulation-repmat(mean(Sat_Simulation,2),1,Num_RML))/sqrt(Num_RML-1);
    
%     if size(Perturb_Flux_Simulation,1) <= size(Perturb_Flux_Simulation,2)
%        verify_ESMDA = verify_ESMDA + Perturb_verify_ESMDA*Perturb_Flux_Simulation'/(Perturb_Flux_Simulation*Perturb_Flux_Simulation'+Iteration_ESMDA*Rq)*(Flux_Perturbation-Flux_Simulation);
%     else
%        xx=(eye(Num_RML)+(Perturb_Flux_Simulation'/(Iteration_ESMDA*Rq))*Perturb_Flux_Simulation)\Perturb_Flux_Simulation'/(Iteration_ESMDA*Rq);
%        verify_ESMDA = verify_ESMDA + Perturb_verify_ESMDA*xx*(Flux_Perturbation-Flux_Simulation);
%     end
%     
%     if size(Perturb_Watercut_Simulation,1) <= size(Perturb_Watercut_Simulation,2)
%        verify_ESMDA = verify_ESMDA + Perturb_verify_ESMDA*Perturb_Watercut_Simulation'/(Perturb_Watercut_Simulation*Perturb_Watercut_Simulation'+Iteration_ESMDA*Rfw)*(Watercut_Perturbation-Watercut_Simulation);
%     else
%        xx=(eye(Num_RML)+(Perturb_Watercut_Simulation'/(Iteration_ESMDA*Rfw))*Perturb_Watercut_Simulation)\Perturb_Watercut_Simulation'/(Iteration_ESMDA*Rfw);
%        verify_ESMDA = verify_ESMDA + Perturb_verify_ESMDA*xx*(Watercut_Perturbation-Watercut_Simulation);
%     end
%     
%     if size(Perturb_BHP_Simulation,1) <= size(Perturb_BHP_Simulation,2)
%        verify_ESMDA = verify_ESMDA + Perturb_verify_ESMDA*Perturb_BHP_Simulation'/(Perturb_BHP_Simulation*Perturb_BHP_Simulation'+Iteration_ESMDA*Rp)*(BHP_Perturbation-BHP_Simulation);
%     else
%        xx=(eye(Num_RML)+(Perturb_BHP_Simulation'/(Iteration_ESMDA*Rp))*Perturb_BHP_Simulation)\Perturb_BHP_Simulation'/(Iteration_ESMDA*Rp);
%        verify_ESMDA = verify_ESMDA + Perturb_verify_ESMDA*xx*(BHP_Perturbation-BHP_Simulation);
%     end
    if Data_type==1
        xx=inv(Rq)-inv(Rq)*Perturb_Flux_Simulation*inv(eye(Num_RML)+Perturb_Flux_Simulation'*inv(Rq)*Perturb_Flux_Simulation)*Perturb_Flux_Simulation'*inv(Rq);
        verify_ESMDA = verify_ESMDA + Perturb_verify_ESMDA*Perturb_Flux_Simulation'*xx*(Flux_Perturbation-Flux_Simulation);

        xx=inv(Rfw)-inv(Rfw)*Perturb_Watercut_Simulation*inv(eye(Num_RML)+Perturb_Watercut_Simulation'*inv(Rfw)*Perturb_Watercut_Simulation)*Perturb_Watercut_Simulation'*inv(Rfw);
        verify_ESMDA = verify_ESMDA + Perturb_verify_ESMDA*Perturb_Watercut_Simulation'*xx*(Watercut_Perturbation-Watercut_Simulation);

        xx=inv(Rp)-inv(Rp)*Perturb_BHP_Simulation*inv(eye(Num_RML)+Perturb_BHP_Simulation'*inv(Rp)*Perturb_BHP_Simulation)*Perturb_BHP_Simulation'*inv(Rp);
        verify_ESMDA = verify_ESMDA + Perturb_verify_ESMDA*Perturb_BHP_Simulation'*xx*(BHP_Perturbation-BHP_Simulation);
        
    elseif Data_type==2
        xx=inv(Rs)-inv(Rs)*Perturb_Sat_Simulation*inv(eye(Num_RML)+Perturb_Sat_Simulation'*inv(Rs)*Perturb_Sat_Simulation)*Perturb_Sat_Simulation'*inv(Rs);
        verify_ESMDA = verify_ESMDA + Perturb_verify_ESMDA*Perturb_Sat_Simulation'*xx*(Sat_Perturbation-Sat_Simulation);    
    else
        
        xx=inv(Rq)-inv(Rq)*Perturb_Flux_Simulation*inv(eye(Num_RML)+Perturb_Flux_Simulation'*inv(Rq)*Perturb_Flux_Simulation)*Perturb_Flux_Simulation'*inv(Rq);
        verify_ESMDA = verify_ESMDA + Perturb_verify_ESMDA*Perturb_Flux_Simulation'*xx*(Flux_Perturbation-Flux_Simulation);

        xx=inv(Rfw)-inv(Rfw)*Perturb_Watercut_Simulation*inv(eye(Num_RML)+Perturb_Watercut_Simulation'*inv(Rfw)*Perturb_Watercut_Simulation)*Perturb_Watercut_Simulation'*inv(Rfw);
        verify_ESMDA = verify_ESMDA + Perturb_verify_ESMDA*Perturb_Watercut_Simulation'*xx*(Watercut_Perturbation-Watercut_Simulation);

        xx=inv(Rp)-inv(Rp)*Perturb_BHP_Simulation*inv(eye(Num_RML)+Perturb_BHP_Simulation'*inv(Rp)*Perturb_BHP_Simulation)*Perturb_BHP_Simulation'*inv(Rp);
        verify_ESMDA = verify_ESMDA + Perturb_verify_ESMDA*Perturb_BHP_Simulation'*xx*(BHP_Perturbation-BHP_Simulation);

        xx=inv(Rs)-inv(Rs)*Perturb_Sat_Simulation*inv(eye(Num_RML)+Perturb_Sat_Simulation'*inv(Rs)*Perturb_Sat_Simulation)*Perturb_Sat_Simulation'*inv(Rs);
        verify_ESMDA = verify_ESMDA + Perturb_verify_ESMDA*Perturb_Sat_Simulation'*xx*(Sat_Perturbation-Sat_Simulation);         
        
    end

end


% %% STEP 1: running ensemble     
% for jj=1:Num_RML
% 
%     verify2=verify_ESMDA(:,jj);
% 
%     for z=1:nz
%         xx=0;
%         for zz=1:z
%             xx=xx+Layer{zz}.D;
%         end
%         verify1=verify2(xx-Layer{zz}.D+1:xx);
%         K=Layer{z}.mPermeablility+Layer{z}.fala*Layer{z}.TMGL*verify1;
%         PORO=0.25 * (exp(K) ./200).^0.1;
%         ROCK=convertFrom(exp(K),milli*darcy);
%         rock{jj}.perm(1+(z-1)*nx*ny:z*nx*ny,1)=ROCK;
%         rock{jj}.poro(1+(z-1)*nx*ny:z*nx*ny,1)=PORO;   
% 
%     end
% 
% end
% 
% num_peturb = Num_RML;
% 
% tic
% Gradient_Calculation_finitedifference;
% toc
% 
% Flux_Simulation=[];
% Watercut_Simulation=[];
% BHP_Simulation=[];
% 
% 
% % ensemble of simulated measurements
% for jj=1:Num_RML
% 
%     % ensemble of simulated measurements
%     Flux11 = FluxFD(1+(jj-1)*nWP:jj*nWP,:);
%     Watercut11 = WatercutFD(1+(jj-1)*nWP:jj*nWP,:);
%     Pres11 = BHPFD(1+(jj-1)*nWI:jj*nWI,:);
% 
%     Flux_Simulation=[Flux_Simulation  reshape(Flux11,size(Flux11,1)*size(Flux11,2),1)];
%     Watercut_Simulation=[Watercut_Simulation  reshape(Watercut11,size(Watercut11,1)*size(Watercut11,2),1)];
%     BHP_Simulation=[BHP_Simulation  reshape(Pres11,size(Pres11,1)*size(Pres11,2),1)];  
% 
% end
% clear Flux11 Watercut11 Pres11
% 
% 
% % ensemble of perturbed measurement
% Flux_Perturbation1 = reshape(Fluxref,size(Fluxref,1)*size(Fluxref,2),1);
% Watercut_Perturbation1 = reshape(Watercutref,size(Watercutref,1)*size(Watercutref,2),1);
% BHP_Perturbation1 = reshape(BHPref,size(BHPref,1)*size(BHPref,2),1);  
% 
% Flux_Perturbation=[];
% Watercut_Perturbation=[];
% BHP_Perturbation=[];
% for jj=1:Num_RML
% 
%     Flux_Perturbation=[Flux_Perturbation  normrnd(Flux_Perturbation1,Flux_Perturbation1*0.05*sqrt(Iteration_ESMDA),size(Fluxref,1)*size(Fluxref,2),1)];
%     Watercut_Perturbation=[Watercut_Perturbation  normrnd(Watercut_Perturbation1, ones(size(Watercutref,1)*size(Watercutref,2),1)*0.05*sqrt(Iteration_ESMDA),size(Watercutref,1)*size(Watercutref,2),1)];
%     BHP_Perturbation=[BHP_Perturbation normrnd(BHP_Perturbation1,BHP_Perturbation1*0.05*sqrt(Iteration_ESMDA),size(BHPref,1)*size(BHPref,2),1)];
% 
% end
% clear Flux_Perturbation1 Watercut_Perturbation1 BHP_Perturbation1
% 
% 
% xx=0;
% for jj=1:Num_RML
% 
%     xx=0;
%     xx=xx+(Flux_Simulation(:,jj)-Flux_Perturbation(:,jj))'*inv(Rq)*(Flux_Simulation(:,jj)-Flux_Perturbation(:,jj));
%     xx=xx+(Watercut_Simulation(:,jj)-Watercut_Perturbation(:,jj))'*inv(Rfw)*(Watercut_Simulation(:,jj)-Watercut_Perturbation(:,jj));
%     xx=xx+(BHP_Simulation(:,jj)-BHP_Perturbation(:,jj))'*inv(Rp)*(BHP_Simulation(:,jj)-BHP_Perturbation(:,jj));
% 
%     RMSE_ensemble_ESMDA(Ite,jj)=xx;
% end
% 
% for jj=1:Num_RML
%     RMSE_K_ESMDA(Ite,jj)=norm(log(convertTo(rock{jj}.perm,milli*darcy))-K_ref,2)/sqrt(length(K_ref));
% end 
