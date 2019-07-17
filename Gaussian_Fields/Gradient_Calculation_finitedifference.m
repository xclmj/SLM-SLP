clear wellDataSGSFD
clear SeismicData1
for ie = 1 : num_peturb
    % well structure
    W = [];
    for iw = 1 : length(wells)
         X = wells{iw}{1}; Y = wells{iw}{2}; Z = 1:G.cartDims(3);
         type = 'rate'; val = wells{iw}{3}; 
         if wells{iw}{3}<0, type = 'bhp'; val = wells{iw}{4}; end
         radius = wells{iw}{5}; comp = wells{iw}{6}; name = wells{iw}{7};
         W = verticalWell(W, G, rock{ie}, X, Y, Z, 'Type', type, ...
             'Val', val, 'Radius', radius, 'Comp_i', comp, 'name', name);
    end
    
    well = W; clear W  
    
    
    if Num_kr==6 && Both_K_Kr==2
        
%         verify=tan(verify_kr(:,ie)*pi/2);
        if ie==num_peturb
           verify_kr=verify_Linearization(Num_K+1:Num_K+Num_kr,1); 
        elseif ie>size(Perturb_Matrix,1)
           verify_kr=verify_Linearization(Num_K+1:Num_K+Num_kr,1)+Perturb_Matrix_kr(:,ie-size(Perturb_Matrix,1));
        else
           verify_kr=verify_Linearization(Num_K+1:Num_K+Num_kr,1); 
        end
        
        verify_kr=verify_kr.*([n_max; n_max; S_max; S_max; kr_max; kr_max]-[n_min; n_min; S_min; S_min; kr_min; kr_min])+[n_min; n_min; S_min; S_min; kr_min; kr_min];
%         verify_kr=(exp(verify_kr).*[n_max; n_max; S_max; S_max; kr_max; kr_max]+[n_min; n_min; S_min; S_min; kr_min; kr_min])./(ones(Num_kr,1)+exp(verify_kr));
        
        fluid = initCoreyFluid('mu' , [ 0.4, 2]*centi*poise     , ...
                               'rho', [1014, 859]*kilogram/meter^3, ...
                               'n'  , [   verify_kr(1),   verify_kr(2)]                 , ...
                               'sr' , [ verify_kr(3), verify_kr(4)]                 , ...
                               'kwm', [  verify_kr(5),   verify_kr(6)]);

    elseif Num_kr==3 && Both_K_Kr==2

%         verify=tan(verify_kr(:,ie)*pi/2);
        if ie==num_peturb
           verify_kr=verify_Linearization(Num_K+1:Num_K+Num_kr,1); 
        elseif ie>size(Perturb_Matrix,1)
           verify_kr=verify_Linearization(Num_K+1:Num_K+Num_kr,1)+Perturb_Matrix_kr(:,ie-size(Perturb_Matrix,1));
        else
           verify_kr=verify_Linearization(Num_K+1:Num_K+Num_kr,1); 
        end
        
%         verify_kr=(exp(verify_kr).*[n_max; n_max; kr_max]+[n_min; n_min; kr_min])./(ones(Num_kr,1)+exp(verify_kr));
        verify_kr=verify_kr.*([n_max; n_max; S_max; S_max; kr_max; kr_max]-[n_min; n_min; S_min; S_min; kr_min; kr_min])+[n_min; n_min; S_min; S_min; kr_min; kr_min];
        
        fluid = initCoreyFluid('mu' , [ 0.4, 2]*centi*poise     , ...
                               'rho', [1014, 859]*kilogram/meter^3, ...
                               'n'  , [   verify_kr(1),   verify_kr(2)]                 , ...
                               'sr' , [ S_mean, S_mean]                 , ...
                               'kwm', [   verify_kr(3),   kr_mean]);     
                           
    else
        
        fluid = initCoreyFluid('mu' , [ 0.4, 2]*centi*poise     , ...
                               'rho', [1014, 859]*kilogram/meter^3, ...
                               'n'  , [ n_mean, n_mean]                 , ...
                               'sr' , [ S_mean, S_mean]                 , ...
                               'kwm', [ kr_mean,  kr_mean]);           
        
    end
    

    % grid transmissibilities
    trans = computeTrans(G, rock{ie});

    % initial state
    rSol = initState(G, well, p0, s0);

    % solve initial pressure
    rSol = incompTPFA(rSol, G, trans, fluid, 'wells', well, ...
                 'LinSolve', linsolve_p);

    %% Time stepping
    t = 0; time = convertTo(t,day); istep = 0; timeSteps=[];
    
    while t < Tend
        
        step = dT; if t + dT > Tend, step = Tend-t; end

        % update saturations
        rSol = implicitTransport(rSol, G, step, rock{ie}, fluid, ...
                         'wells', well, 'LinSolve', linsolve_t);  
        % Check for inconsistent saturations
        assert(max(rSol.s(:,1)) < 1+eps && min(rSol.s(:,1)) > -eps);
        % Update solution of pressure equation
        rSol = incompTPFA(rSol, G, trans, fluid, 'wells', well, ...
                 'LinSolve', linsolve_p);

        t = t + step; time = convertTo(t,day); istep = istep + 1;
        
        % store well rates and pressures
         if ~exist('wellDataSGSFD','var')
             for i = 1 : length(well)
                 wellDataSGSFD(i).oil  = [];
                 wellDataSGSFD(i).bhp  = [];
                 wellDataSGSFD(i).flx  = [];
                 wellDataSGSFD(i).ffl  = [];
             end
         end

         for i = 1 : length(rSol.wellSol)
             wellDataSGSFD(i).bhp(ie,istep) = rSol.wellSol(i).pressure;
             rate = 0;
             for k=1:length(rSol.wellSol(i).flux)
                 rate = rate + rSol.wellSol(i).flux(k);
             end
             wellDataSGSFD(i).flx(ie,istep) = abs(convertTo(-rate, meter^3/day));

             % compute fractional flow assuming horizontal flow and no 
             % capillary pressure: fw = qw/qt = 1/(1+(muw/krw)*(kro/muo))
             s = mean(rSol.s(well(i).cells));
             kr = fluid.relperm(s); % [krw kro]
             mu = fluid.properties(s); % [muw muo]
             mo = bsxfun(@rdivide, kr, mu); % [krw/muw kro/muo]
             fw = 1/(1 + mo(2)/mo(1));
             if isnan(fw), fw = 0; end
             wellDataSGSFD(i).ffl(ie,istep) = fw;
             wellDataSGSFD(i).oil(ie,istep) = wellDataSGSFD(i).flx(ie,istep)*(1-fw);
         end

         clear rate
         
         SeismicData1(G.cells.indexMap,istep,ie)=rSol.s(:,1);

    end

end

FluxFD=[];
WatercutFD=[];
BHPFD=[];
SatFD=[];

Fluxm=[];
Watercutm=[];
BHPm=[];
ddd=0;

Flux1=[];
Watercut1=[];
BHP1=[];
Sat11=[];

for ii=1:num_peturb
    
    for jj=1:nt

        if mod(jj,mm_welldata)==0
    
           for j=1:1:nWI
               BHP1=[BHP1 ;  wellDataSGSFD(j).bhp(ii,jj)];
           end   
           for j=1:1:nWP
               Watercut1=[Watercut1 ;  wellDataSGSFD(j+nWI).ffl(ii,jj)];
               Flux1=[Flux1 ;  wellDataSGSFD(j+nWI).flx(ii,jj)];
           end   
        end
        
        Fluxm=[Fluxm Flux1];
        Watercutm=[Watercutm Watercut1];
        BHPm=[BHPm BHP1];
        
        Flux1=[];
        Watercut1=[];
        BHP1=[];
        
        if mod(jj,mm_seismicdata)==0
           Sat11=[Sat11 SeismicData1(:,jj,ii)];      %% I select the standard derivation is 10% of real observation 
        end
        

    end
    
    FluxFD=[FluxFD ; Fluxm];
    WatercutFD=[WatercutFD ;  Watercutm];
    BHPFD=[BHPFD ;  BHPm];
    SatFD=[SatFD ; Sat11];
    
    Fluxm=[];
    Watercutm=[];
    BHPm=[];
    Sat11=[];
    
end
clear Flux1 Watercut1 BHP1 Fluxm Watercutm BHPm Sat11