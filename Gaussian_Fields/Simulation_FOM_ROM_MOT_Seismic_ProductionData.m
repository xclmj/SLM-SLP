clear wellDataSGSFD
clear SeismicData1

if FOM_ROM==1
    
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

             SeismicData1(G.cells.indexMap,istep)=rSol.s(:,1);

        end

    end
    
    SatFD=[];

    Sat11=[];

    for ii=1:num_peturb

        for jj=1:nt

            if mod(jj,mm_seismicdata)==0
               Sat11=[Sat11 SeismicData1(:,jj)];      %% I select the standard derivation is 10% of real observation 
            end


        end

        SatFD=[SatFD ; Sat11];

        Sat11=[];

    end
    clear Flux1 Watercut1 BHP1 Fluxm Watercutm BHPm Sat11
    
elseif FOM_ROM==2
    
    SatFD=[];
   
    for ISTEP=1:ddd_seismicdata

        SatFD(:,ISTEP)=Seismic_ROLM(:,ISTEP)+Gra_Seismic2(:,:,ISTEP)*(verify-verify_Linearization);

    end
    
    FluxFD=[];
    WatercutFD=[];
    BHPFD=[];
    
    for ISTEP=1:ddd_welldata

        FluxFD(:,ISTEP)=Flux_ROLM(:,ISTEP)+Gra_Flux2(:,:,ISTEP)*(verify-verify_Linearization);
        WatercutFD(:,ISTEP)=Watercut_ROLM(:,ISTEP)+Gra_Wat2(:,:,ISTEP)*(verify-verify_Linearization);
        BHPFD(:,ISTEP)=BHP_ROLM(:,ISTEP)+Gra_BHP2(:,:,ISTEP)*(verify-verify_Linearization);
        
    end


end