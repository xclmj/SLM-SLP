%% step 1: generate pertubated measurement with reference
mm_welldata=2;                % timesteps for well data
mm_seismicdata=25;            % timesteps for seismic data

%
K_ref=[];            % reference log-permeability field
K=[];
for z=1:nz
    
    verify=(rand(size(Layer{z}.fala,2),1)-0.5*ones(size(Layer{z}.fala,2),1))*1.5;
    K=Layer{z}.mPermeablility+Layer{z}.fala*verify;
    verify_k_ref=verify;
    verify_G=(rand(nx*ny,1)-0.5*ones(nx*ny,1))*1.5;
    verify_G(1:length(verify))=verify;

    PORO=0.25 * (exp(K) ./200).^0.1;
    ROCK=convertFrom(exp(K),milli*darcy);
    rock{1}.perm(1+(z-1)*nx*ny:z*nx*ny,1)=ROCK;
    rock{1}.poro(1+(z-1)*nx*ny:z*nx*ny,1)=PORO;     

end

if Num_kr==6  % reference relative permeability parameters
   verify_kr_ref(1)= n_mean; verify_kr_ref(2)= n_mean;
   verify_kr_ref(3)= S_mean; verify_kr_ref(4)= S_mean;       
   verify_kr_ref(5)= kr_mean; verify_kr_ref(6)= kr_mean;   
else
   verify_kr_ref(1)= n_mean; verify_kr_ref(2)= n_mean;
   verify_kr_ref(3)= kr_mean;       
end
%  if Num_kr==6  % reference relative permeability parameters
%    verify_kr_ref(1)= n_max; verify_kr_ref(2)= n_max;
%    verify_kr_ref(3)= S_max; verify_kr_ref(4)= S_max;       
%    verify_kr_ref(5)= kr_max; verify_kr_ref(6)= kr_max;   
% else
%    verify_kr_ref(1)= n_max; verify_kr_ref(2)= n_max;
%    verify_kr_ref(3)= kr_max;       
% end

K_ref=log(convertTo(rock{1}.perm,milli*darcy));
ValidationSolution;

s = linspace(0, 1, 1001).'; kr = fluid.relperm(s);
plot(s, kr), legend('kr_1', 'kr_2'), hold on

%%%%%%----generating observation perturbation for MAP Problem
%1: assimilation well data
Rq=[];            % covariance for production flux
Rfw=[];           % covariance for water cut
Rp=[];            % covariance for bottom-hole pressure at injection wells

%2:assimilation seismic data (using grid saturation)
Rs=[];            % covariance for saturation at every grid

Flux1=[];
Watercut1=[];
Pres1=[];
Sat1=[];

% ensemble of measurement
Flux=[];
Watercut=[];
BHP=[];
Sat=[];

ddd_welldata=0;               
ddd_seismicdata=0;           
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
    
    if mod(j,mm_seismicdata)==0
       ddd_seismicdata=ddd_seismicdata+1;
       Rs(1+(ddd_seismicdata-1)*nx*ny*nz:ddd_seismicdata*nx*ny*nz,1+(ddd_seismicdata-1)*nx*ny*nz:ddd_seismicdata*nx*ny*nz)=diag((SeismicData(:,j)*0.1).^2);      %% I select the standard derivation is 10% of real observation 
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
    
    if mod(jj,mm_seismicdata)==0
       Sat=[Sat normrnd(SeismicData(:,jj),SeismicData(:,jj)*0.1,nx*ny*nz,1)];      %% I select the standard derivation is 10% of real observation 
    end
    
    
end
clear Flux1 Watercut1 Pres1 Sat1

%%%%%%----generating observation perturbation for RML Problem
Flux1=[];
Watercut1=[];
Pres1=[];

Flux11=[];
Watercut11=[];
Pres11=[];
Sat11=[];

% ensemble of measurement
Flux_RML=[];
Watercut_RML=[];
BHP_RML=[];
Sat_RML=[];

Num_RML=500;
for kk=1:Num_RML
    
    for jj=1:nt

        if mod(jj,mm_welldata)==0     
           for j=1:1:nWI
               Pres1=[Pres1 ; normrnd(wellDataSGSo(j).bhp(1,jj),wellDataSGSo(j).bhp(1,jj)*0.1)];
           end  

           for j=1:1:length(rSol.wellSol)-nWI
               Flux1=[Flux1 ; normrnd(wellDataSGSo(j+nWI).flx(1,jj),wellDataSGSo(j+nWI).flx(1,jj)*0.1)];
               Watercut1=[Watercut1 ; normrnd(wellDataSGSo(j+nWI).ffl(1,jj),0.05)];
           end       

        end

        Flux11=[Flux11 Flux1];
        Watercut11=[Watercut11 Watercut1];
        Pres11=[Pres11 Pres1];

        Flux1=[];
        Watercut1=[];
        Pres1=[];
        
        if mod(jj,mm_seismicdata)==0
           Sat11=[Sat11 normrnd(SeismicData(:,jj),SeismicData(:,jj)*0.1,nx*ny*nz,1)];      %% I select the standard derivation is 10% of real observation 
        end

    end
    
    Flux_RML=[Flux_RML ; Flux11];
    Watercut_RML=[Watercut_RML ; Watercut11];
    BHP_RML=[BHP_RML ; Pres11];  
    Sat_RML=[Sat_RML ; Sat11];

    Flux11=[];
    Watercut11=[];
    Pres11=[];  
    Sat11=[];

end
clear Flux1 Watercut1 Pres1 Sat11 Flux11 Watercut11 Pres11

% initial ensemble parameters
Num_K=0;
for z=1:nz
    Num_K=Num_K+sum(Layer{z}.falaK);
end

for ii=1:Num_RML
    
    verify_RML_Init(1:Num_K,ii)=(rand(Num_K,1)-0.5*ones(Num_K,1))*2;
    
end

%% real simulated measurement at reference permeability
Fluxref=[];
Watercutref=[];
BHPref=[];
Satref=[];

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
    
    if mod(jj,mm_seismicdata)==0
       Satref=[Satref SeismicData(:,jj)];      %% I select the standard derivation is 10% of real observation 
    end
    
end
clear Flux1 Watercut1 Pres1
