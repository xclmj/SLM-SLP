%% Smolyak Sparse Grid Solution
 %Tsim=5*year();
clear wellDataSGSPC_FD
%     rock{1}.perm=K_uptFD;
%     rock{2}.perm=K_ref;
%     rock{3}.perm=mPermeablility;
for ie = 1 : 1   %size(update_ref,2)+1
    
%     if ie==1
%        rock{ie}.perm=mPermeablility;
%     else
%        rock{ie}.perm=mPermeablility+fala*update_ref(:,ie-1);
%     end
    
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

%     VSsaturation(:,1,ie) = rSol.s(:,2);
%     VSpressure(:,1,ie) = rSol.pressure;         
    %% Time stepping
    t = 0; time = convertTo(t,day); istep = 0; timeSteps=[];
    
    tic
    while t < Tsim
        
        step = dT; if t + dT > Tsim, step = Tsim-t; end

        % update saturations
        rSol = implicitTransport(rSol, G, step, rock{ie}, fluid, ...
                         'wells', well, 'LinSolve', linsolve_t);  
        % Check for inconsistent saturations
        assert(max(rSol.s(:,1)) < 1+eps && min(rSol.s(:,1)) > -eps);
        % Update solution of pressure equation
        rSol = incompTPFA(rSol, G, trans, fluid, 'wells', well, ...
                 'LinSolve', linsolve_p);

        t = t + step; time = convertTo(t,day); istep = istep + 1;

        storeWellDataSGSPC_FD;
        % store snapshots of saturation and pressure
        
        if ie == 1; timeSteps = [timeSteps; time]; end

    end
    toc
end