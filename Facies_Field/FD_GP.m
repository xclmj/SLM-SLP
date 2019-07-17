%% Gradient-based History Matching by finite difference
Jflux=0;
Jwatercut=0;
Jpres=0;
JtotalFD_G=[];

% global cost function at the reference field
for tt=1:ddd_welldata
    tttt1=(1+(tt-1)*nWI):(nWI+(tt-1)*nWI);
    Jpres=Jpres+(BHPref(:,tt)-BHP(:,tt))'*inv(Rp(tttt1,tttt1))*(BHPref(:,tt)-BHP(:,tt))/2;
    tttt2=(1+(tt-1)*nWP):(nWP+(tt-1)*nWP);  
    Jwatercut=Jwatercut+(Watercutref(:,tt)-Watercut(:,tt))'*inv(Rfw(tttt2,tttt2))*(Watercutref(:,tt)-Watercut(:,tt))/2; 
    Jflux=Jflux+(Fluxref(:,tt)-Flux(:,tt))'*inv(Rq(tttt2,tttt2))*(Fluxref(:,tt)-Flux(:,tt))/2;
end    
Jtotal=Jflux+Jwatercut+Jpres;
JtotalFD_G=[JtotalFD_G Jtotal];

EnK_num=5;
Domain_Numx=3;
Domain_Numy=3;     
Initilization_DD1;
parameterization_new1;
Num_K=0;

Num_K=0;
for z=1:nz
    Num_K=Num_K+size(Layer{z}.fala,2);
end
num_peturb=Num_K+1;

verify=zeros(Num_K,1);
%verify=(rand(dd,1)-0.5*ones(dd,1));
erro1=1e-4;  % iteration criteration
erro2=1e-3;
IteNum=100;
sigma1=1;
sigma2=1;


Stop_Condition=0;

for JJ=1:IteNum
    tic
    for jj=1:num_peturb

        if jj<num_peturb
           XX=zeros(Num_K,1); XX(jj)=0.01;
           verify2=verify+XX;
        else
           verify2=verify;
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
    
    JtotalFD_G=[JtotalFD_G JtotalFDm_G(end)];
    JtotalFDm_G(end)
     
    % direct computation of global gradient
    b=[];
    A=[];
    Gra_G=[];
    for i=1:num_peturb-1
        Gra_G(i,1)=(JtotalFDm_G(i)-JtotalFDm_G(num_peturb))/0.01;
    end
    
    if JJ>1   

        sigma1=abs(JtotalFD_G(end)-JtotalFD_G(end-1))/max([1 JtotalFD_G(end)]);
%             sigma2=norm(verify-verifym,2)/max([1 norm(verify,2)]);
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
    
    EE=Gra_G/max(abs(Gra_G));
    
    if JJ==1
       ll=0.2;
    end

   if  JJ>1 && JtotalFD_G(end)>JtotalFD_G(end-1)
        ll=ll/2;
        verify=verifym-ll*FF;
    else 

         verifym=verify;
         verify=verifym-ll*EE;               
    end

    FF=EE;
    verify_upt=verify;
    AAA=['Inner-loop Iteration Step:' num2str(JJ)];
    disp(AAA)

end

for z=1:nz
    xx=0;
    for zz=1:z
        xx=xx+size(Layer{zz}.fala,2);
    end
    verify1=verify(xx-size(Layer{zz}.fala,2)+1:xx);
    K_upt=Layer{z}.mPermeablility+Layer{z}.fala*verify1;
    [K_upt,S]= OPCA_TwoFacies(Layer{z}.fala, verify1, Layer{z}.mPermeablility, gammaF);
    K_uptFD_G=K_upt;
end 

figure(1)
plot(1:length(JtotalFD_G)-1,JtotalFD_L(2:length(JtotalFD_G)),'b')
hold on
plot(1:length(JtotalFD_G)-1,JtotalFD_L(1)*ones(length(JtotalFD_G)-1,1),'r','LineWidth',2);

xlabel('Iteration Step')
ylabel('Cost Function')
