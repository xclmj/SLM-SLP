%% Method1: without graph coloring

if nz==1 && DP_IP==1
    % Perturbation 1: as Fig.2 in SPE 130473 "Development of a Data Partition Technique for
    % Gradient-Based Optimization Methods in History Matching"
    delta=0.01;
    n1=0; n2=0; n3=0; n4=0; n5=0;
    clear Perturb_Index Perturb_Matrix

    for j=1:Domain_Num

        if rem(Domain_index(j,1),3)==1
           n1=n1+1;
           Perturb_Index{1}.ID(n1)=j;
        elseif rem(Domain_index(j,1),3)==2
           n2=n2+1;
           Perturb_Index{2}.ID(n2)=j;        
        elseif rem(Domain_index(j,1),3)==0
           n5=n5+1;
           Perturb_Index{5}.ID(n5)=j;  
        end

        if rem(Domain_index(j,2),3)==1
           n3=n3+1;
           Perturb_Index{3}.ID(n3)=j;
        elseif rem(Domain_index(j,2),3)==2
           n4=n4+1;
           Perturb_Index{4}.ID(n4)=j;        
        elseif rem(Domain_index(j,2),3)==0
           n5=n5+1;
           Perturb_Index{5}.ID(n5)=j;  
        end    

    end

    Perturb_Matrix=zeros(5*EnK_num,EnK_num*Domain_Num);
    for i=1:5
        tt1=1+(i-1)*EnK_num:(EnK_num+(i-1)*EnK_num);        
        for j=1:length(Perturb_Index{i}.ID)
            tt2=1+(Perturb_Index{i}.ID(j)-1)*EnK_num:(EnK_num+(Perturb_Index{i}.ID(j)-1)*EnK_num);
            Perturb_Matrix(tt1,tt2)=delta*eye(EnK_num,EnK_num);
        end  
    end

    % % Perturbation 2: Vandermonde matrix
    % delta=[];
    % for i=1:Domain_Num
    %     delta=[delta 0.1+0.01*(i-1)];
    % end
    % Perturb_Matrix=zeros(5*EnK_num,EnK_num*Domain_Num);
    % for i=1:5
    %     tt1=1+(i-1)*EnK_num:(EnK_num+(i-1)*EnK_num);
    %     for j=1:Domain_Num
    %         tt2=1+(j-1)*EnK_num:(EnK_num+(j-1)*EnK_num);
    %         Perturb_Matrix(tt1,tt2)=delta(j)^(i-1)*eye(EnK_num,EnK_num);
    %     end
    % end
    % for i=1:EnK_num*Domain_Num
    %     Perturb_Matrix(:,i)=Perturb_Matrix(:,i)/sqrt(Perturb_Matrix(:,i)'*Perturb_Matrix(:,i));
    % end
    % 
    % % Perturbation 3: Chebyshev matrix
    % Perturb_Matrix=zeros(5*EnK_num,EnK_num*Domain_Num);
% for i=1:5
%     tt1=1+(i-1)*EnK_num:(EnK_num+(i-1)*EnK_num);
%     for j=1:Domain_Num
%         tt2=1+(j-1)*EnK_num:(EnK_num+(j-1)*EnK_num);
%         Perturb_Matrix(tt1,tt2)=cos((2*j-1)*i*pi/2/Domain_Num)*eye(EnK_num,EnK_num);
%     end  
% end
% for i=1:EnK_num*Domain_Num
%     Perturb_Matrix(:,i)=Perturb_Matrix(:,i)/sqrt(Perturb_Matrix(:,i)'*Perturb_Matrix(:,i));
% end

%% Method2: with graph coloring
% referecn to the SPE 140811 "Study of Perturbation Designs for the Partially 
% Separable Objective Function With Gradient-Based Optimizations in History Matching"
% %% Method3: independent perturbation for each subdomain

elseif nz==2 && Vandermonde_Chebyshev==1 && DP_IP==1
    
    clear Perturb_Index Perturb_Matrix

    % % Perturbation 2: Vandermonde matrix
    delta=[];
    for i=1:Domain_Num*nz
        delta=[delta 0.1+0.01*(i-1)];
    end
    Perturb_Matrix=zeros(6*EnK_num,EnK_num*Domain_Num*nz);
    for i=1:6
        tt1=1+(i-1)*EnK_num:(EnK_num+(i-1)*EnK_num);
        for j=1:Domain_Num*nz
            tt2=1+(j-1)*EnK_num:(EnK_num+(j-1)*EnK_num);
            Perturb_Matrix(tt1,tt2)=delta(j)^(i-1)*eye(EnK_num,EnK_num);
        end
    end
    for i=1:EnK_num*Domain_Num*nz
        Perturb_Matrix(:,i)=Perturb_Matrix(:,i)/sqrt(Perturb_Matrix(:,i)'*Perturb_Matrix(:,i));
    end
    Perturb_Matrix=Perturb_Matrix/10;
    
elseif nz==2 && Vandermonde_Chebyshev==2 && DP_IP==1
    
    clear Perturb_Index Perturb_Matrix
    % % Perturbation 3: Chebyshev matrix
    Perturb_Matrix=zeros(6*EnK_num,EnK_num*Domain_Num*nz);
    for i=1:6
        tt1=1+(i-1)*EnK_num:(EnK_num+(i-1)*EnK_num);
        for j=1:Domain_Num*nz
            tt2=1+(j-1)*EnK_num:(EnK_num+(j-1)*EnK_num);
            Perturb_Matrix(tt1,tt2)=cos((2*j-1)*i*pi/2/Domain_Num)*eye(EnK_num,EnK_num);
        end  
    end
    for i=1:EnK_num*Domain_Num*nz
        Perturb_Matrix(:,i)=Perturb_Matrix(:,i)/sqrt(Perturb_Matrix(:,i)'*Perturb_Matrix(:,i));
    end
    Perturb_Matrix=Perturb_Matrix/10;

elseif nz>2 && Vandermonde_Chebyshev==1 && DP_IP==1
    
    clear Perturb_Index Perturb_Matrix

    % % Perturbation 2: Vandermonde matrix
    delta=[];
    for i=1:Domain_Num*nz
        delta=[delta 0.1+0.01*(i-1)];
    end
    Perturb_Matrix=zeros(7*EnK_num,EnK_num*Domain_Num*nz);
    for i=1:7
        tt1=1+(i-1)*EnK_num:(EnK_num+(i-1)*EnK_num);
        for j=1:Domain_Num*nz
            tt2=1+(j-1)*EnK_num:(EnK_num+(j-1)*EnK_num);
            Perturb_Matrix(tt1,tt2)=delta(j)^(i-1)*eye(EnK_num,EnK_num);
        end
    end
    for i=1:EnK_num*Domain_Num*nz
        Perturb_Matrix(:,i)=Perturb_Matrix(:,i)/sqrt(Perturb_Matrix(:,i)'*Perturb_Matrix(:,i));
    end
    Perturb_Matrix=Perturb_Matrix/10;
    
elseif nz>2 && Vandermonde_Chebyshev==2 && DP_IP==1
    
    clear Perturb_Index Perturb_Matrix
    % % Perturbation 3: Chebyshev matrix
    Perturb_Matrix=zeros(7*EnK_num,EnK_num*Domain_Num*nz);
    for i=1:7
        tt1=1+(i-1)*EnK_num:(EnK_num+(i-1)*EnK_num);
        for j=1:Domain_Num*nz
            tt2=1+(j-1)*EnK_num:(EnK_num+(j-1)*EnK_num);
            Perturb_Matrix(tt1,tt2)=cos((2*j-1)*i*pi/2/Domain_Num)*eye(EnK_num,EnK_num);
        end  
    end
    for i=1:EnK_num*Domain_Num*nz
        Perturb_Matrix(:,i)=Perturb_Matrix(:,i)/sqrt(Perturb_Matrix(:,i)'*Perturb_Matrix(:,i));
    end
    Perturb_Matrix=Perturb_Matrix/10;
    
elseif DP_IP==2
    
    delta=0.01;
    clear Perturb_Index Perturb_Matrix

    Perturb_Matrix=zeros(nz*EnK_num,EnK_num*Domain_Num*nz);
    for i=1:3
        tt1=1+(i-1)*EnK_num:(EnK_num+(i-1)*EnK_num);        
        for j=1:length(Perturb_Index{i}.ID)
            tt2=1+(Perturb_Index{i}.ID(j)-1)*EnK_num:(EnK_num+(Perturb_Index{i}.ID(j)-1)*EnK_num);
            Perturb_Matrix(tt1,tt2)=delta*eye(EnK_num,EnK_num);
        end  
    end  
    
 
end








