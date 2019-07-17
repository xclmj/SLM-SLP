
disp('Measurement Generation')
%% step 1: generate pertubated measurement with reference
mm_welldata=2;                % timesteps for well data
mm_seismicdata=25;            % timesteps for seismic data

K_ref=[];            % reference field
ROCK=[];
PORO=[];
for z=1:nz   
    verify=(rand(size(Layer{z}.fala,2),1)-0.5*ones(size(Layer{z}.fala,2),1))*2;
    K=Layer{z}.mPermeablility+Layer{z}.fala*verify;
    K_ref_LPCA=K;  
    [K,S]= OPCA_TwoFacies(Layer{z}.fala, verify, Layer{z}.mPermeablility, gammaF);
    K1=K_mud*exp(log(K_sand/K_mud)*K); 
    PORO=[PORO ; 0.25 * (K1 ./200).^0.1];
    ROCK=[ROCK ; convertFrom(K1,milli*darcy)];       
end

rock{1}.perm(G.cells.indexMap,1)=ROCK(G.cells.indexMap);
rock{1}.poro(G.cells.indexMap,1)=PORO(G.cells.indexMap); 
    
K_ref=K;
ValidationSolution;

%%%%%%----generating observation perturbation for MAP Problem
%1: assimilation well data
Rq=[];            % covariance for production flux
Rfw=[];           % covariance for water cut
Rp=[];            % covariance for bottom-hole pressure at injection wells

Flux1=[];
Watercut1=[];
Pres1=[];

% ensemble of measurement
Flux=[];
Watercut=[];
BHP=[];

ddd_welldata=0;                        
for j=1:nt
    
    if mod(j,mm_welldata)==0
       ddd_welldata=ddd_welldata+1;
       for i=1:nWI
           Rp(i+(ddd_welldata-1)*nWI,i+(ddd_welldata-1)*nWI)=(wellDataSGSo(i).bhp(1,j)*0.05)^2;      %% I select the standard derivation is 10% of real observation 
       end
 
       for i=1:nWP
           Rfw(i+(ddd_welldata-1)*nWP,i+(ddd_welldata-1)*nWP)=(0.05)^2;      %% I select the standard derivation is 10% of real observation 
           Rq(i+(ddd_welldata-1)*nWP,i+(ddd_welldata-1)*nWP)=(wellDataSGSo(i+nWI).flx(1,j)*0.05)^2;
       end

    end
    
end


for jj=1:nt
    
    if mod(jj,mm_welldata)==0    
       for j=1:1:nWI
           Pres1=[Pres1 ; normrnd(wellDataSGSo(j).bhp(1,jj),wellDataSGSo(j).bhp(1,jj)*0.05)];
       end  
       
       for j=1:1:length(rSol.wellSol)-nWI
           Flux1=[Flux1 ; normrnd(wellDataSGSo(j+nWI).flx(1,jj),wellDataSGSo(j+nWI).flx(1,jj)*0.05)];
           Watercut1=[Watercut1 ; normrnd(wellDataSGSo(j+nWI).ffl(1,jj),0.05)];
       end       

    end
    
    Flux=[Flux Flux1];
    Watercut=[Watercut Watercut1];
    BHP=[BHP Pres1];
    
    Flux1=[];
    Watercut1=[];
    Pres1=[];

end
clear Flux1 Watercut1 Pres1

%%%%%%----generating observation perturbation for RML Problem
Flux1=[];
Watercut1=[];
Pres1=[];

Flux11=[];
Watercut11=[];
Pres11=[];

% ensemble of measurement
Flux_RML=[];
Watercut_RML=[];
BHP_RML=[];

Num_RML=30;
for kk=1:Num_RML
    
    for jj=1:nt

        if mod(jj,mm_welldata)==0     
           for j=1:1:nWI
               Pres1=[Pres1 ; normrnd(wellDataSGSo(j).bhp(1,jj),wellDataSGSo(j).bhp(1,jj)*0.05)];
           end  

           for j=1:1:length(rSol.wellSol)-nWI
               Flux1=[Flux1 ; normrnd(wellDataSGSo(j+nWI).flx(1,jj),wellDataSGSo(j+nWI).flx(1,jj)*0.05)];
               Watercut1=[Watercut1 ; normrnd(wellDataSGSo(j+nWI).ffl(1,jj),0.05)];
           end       

        end

        Flux11=[Flux11 Flux1];
        Watercut11=[Watercut11 Watercut1];
        Pres11=[Pres11 Pres1];

        Flux1=[];
        Watercut1=[];
        Pres1=[];

    end
    
    Flux_RML=[Flux_RML ; Flux11];
    Watercut_RML=[Watercut_RML ; Watercut11];
    BHP_RML=[BHP_RML ; Pres11];  

    Flux11=[];
    Watercut11=[];
    Pres11=[];  

end
clear Flux1 Watercut1 Pres1 Flux11 Watercut11 Pres11


%% real simulated measurement at reference permeability
Fluxref=[];
Watercutref=[];
BHPref=[];

Flux1=[];
Watercut1=[];
Pres1=[];

for jj=1:nt
    
    if mod(jj,mm_welldata)==0       
       for j=1:1:nWI
           Pres1=[Pres1 ; wellDataSGSo(j).bhp(1,jj)];
       end  
       
       for j=1:1:length(rSol.wellSol)-nWI
           Flux1=[Flux1 ; wellDataSGSo(j+nWI).flx(1,jj)];
           Watercut1=[Watercut1 ; wellDataSGSo(j+nWI).ffl(1,jj)];
       end       

    end
    
    Fluxref=[Fluxref Flux1];
    Watercutref=[Watercutref Watercut1];
    BHPref=[BHPref Pres1];
    
    Flux1=[];
    Watercut1=[];
    Pres1=[];
    
end
clear Flux1 Watercut1 Pres1

%% Or you can load our results used in the paper.
% load Case1.mat
