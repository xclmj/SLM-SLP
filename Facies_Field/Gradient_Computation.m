% computing gradient

if OPtimization_Algorithm==3
    
    EE=verify;
    HH=eye(Num_K,Num_K);

    for tt=1:ddd_welldata

        tttt3=(1+(tt-1)*nWI):(nWI+(tt-1)*nWI);
        tttt4=(1+(tt-1)*nWP):(nWP+(tt-1)*nWP);

        EE=EE+(Gra_BHP(:,:,tt))'*inv(Rp(tttt3,tttt3))*(BHPFD(:,tt)-BHP(:,tt));
        EE=EE+(Gra_Flux(:,:,tt))'*inv(Rq(tttt4,tttt4))*(FluxFD(:,tt)-Flux(:,tt));
        EE=EE+(Gra_Wat(:,:,tt))'*inv(Rfw(tttt4,tttt4))*(WatercutFD(:,tt)-Watercut(:,tt));
        
        HH=HH+(Gra_Flux(:,:,tt))'*inv(Rq(tttt4,tttt4))*Gra_Flux(:,:,tt);
        HH=HH+(Gra_Wat(:,:,tt))'*inv(Rfw(tttt4,tttt4))*(Gra_Wat(:,:,tt));
        HH=HH+(Gra_BHP(:,:,tt))'*inv(Rp(tttt4,tttt4))*(Gra_BHP(:,:,tt));

    end
    
    EE=HH\EE;
    
    
elseif OPtimization_Algorithm==2
    
    if JJ==1
        
        EE=verify;
        HH=eye(Num_K,Num_K);

        for tt=1:ddd_welldata

            tttt3=(1+(tt-1)*nWI):(nWI+(tt-1)*nWI);
            tttt4=(1+(tt-1)*nWP):(nWP+(tt-1)*nWP);

            EE=EE+(Gra_BHP(:,:,tt))'*inv(Rp(tttt3,tttt3))*(BHPFD(:,tt)-BHP(:,tt));
            EE=EE+(Gra_Flux(:,:,tt))'*inv(Rq(tttt4,tttt4))*(FluxFD(:,tt)-Flux(:,tt));
            EE=EE+(Gra_Wat(:,:,tt))'*inv(Rfw(tttt4,tttt4))*(WatercutFD(:,tt)-Watercut(:,tt));

        end
        GG=EE;
        
    else
        
        SS=verify-verifym;
        
        EE=verify;

        for tt=1:ddd_welldata

            tttt3=(1+(tt-1)*nWI):(nWI+(tt-1)*nWI);
            tttt4=(1+(tt-1)*nWP):(nWP+(tt-1)*nWP);

            EE=EE+(Gra_BHP(:,:,tt))'*inv(Rp(tttt3,tttt3))*(BHPFD(:,tt)-BHP(:,tt));
            EE=EE+(Gra_Flux(:,:,tt))'*inv(Rq(tttt4,tttt4))*(FluxFD(:,tt)-Flux(:,tt));
            EE=EE+(Gra_Wat(:,:,tt))'*inv(Rfw(tttt4,tttt4))*(WatercutFD(:,tt)-Watercut(:,tt));

        end
        GG=EE;
        
        YY=EE-WW;
        YY1=SS'*YY+YY'*HH*YY;
        YY1=YY1/(SS'*YY)/(SS'*YY);
        YY1=YY1*SS*SS';
        XX=HH*YY*SS'+SS*(YY)'*HH;
        XX=XX/(SS'*YY);
        HH=HH+YY1-XX;
       
    end
    
    EE=HH*EE;
    
    WW=GG;
    
    
elseif OPtimization_Algorithm==1
    
    EE=verify;

    for tt=1:ddd_welldata

        tttt3=(1+(tt-1)*nWI):(nWI+(tt-1)*nWI);
        tttt4=(1+(tt-1)*nWP):(nWP+(tt-1)*nWP);

        EE=EE+(Gra_BHP(:,:,tt))'*inv(Rp(tttt3,tttt3))*(BHPFD(:,tt)-BHP(:,tt));
        EE=EE+(Gra_Flux(:,:,tt))'*inv(Rq(tttt4,tttt4))*(FluxFD(:,tt)-Flux(:,tt));
        EE=EE+(Gra_Wat(:,:,tt))'*inv(Rfw(tttt4,tttt4))*(WatercutFD(:,tt)-Watercut(:,tt));

    end  
    
    HH=eye(Num_K,Num_K);
    
    EE=EE/max(abs(EE));
    
end
    
    