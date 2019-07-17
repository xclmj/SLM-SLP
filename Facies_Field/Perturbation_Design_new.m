%% Method1: without graph coloring

if Vandermonde_Chebyshev>2
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

    XX=zeros(5*EnK_num*nz,EnK_num*Domain_Num*nz);
    for i=1:5
        tt1=1+(i-1)*EnK_num*nz:(EnK_num*nz+(i-1)*EnK_num*nz);        
        for j=1:length(Perturb_Index{i}.ID)
            tt2=1+(Perturb_Index{i}.ID(j)-1)*EnK_num*nz:(EnK_num*nz+(Perturb_Index{i}.ID(j)-1)*EnK_num*nz);
            XX(tt1,tt2)=delta*eye(EnK_num*nz,EnK_num*nz);
        end  
    end
    
    Perturb_Matrix=zeros(5*EnK_num*nz,EnK_num*Domain_Num*nz);
    for z=1:nz
        for j=1:Domain_Num
            tt1=1+(j-1)*EnK_num+(z-1)*EnK_num*Domain_Num:j*EnK_num+(z-1)*EnK_num*Domain_Num;
            tt2=nz*EnK_num*(j-1)+1+(z-1)*EnK_num:nz*EnK_num*(j-1)+z*EnK_num;
             Perturb_Matrix(1:5*EnK_num*nz,tt1)=XX(1:5*EnK_num*nz,tt2);
        end
    end


elseif Vandermonde_Chebyshev==1 
    
    clear Perturb_Index Perturb_Matrix

    % % Perturbation 2: Vandermonde matrix
    delta=[];
    for i=1:Domain_Num*nz
        delta=[delta 0.1+0.01*(i-1)];
    end
    Perturb_Matrix=zeros(5*EnK_num*nz,EnK_num*Domain_Num*nz);
    for i=1:5
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
    
elseif Vandermonde_Chebyshev==2 
    
    clear Perturb_Index Perturb_Matrix
    % % Perturbation 3: Chebyshev matrix
    Perturb_Matrix=zeros(5*EnK_num*nz,EnK_num*Domain_Num*nz);
    for i=1:5*nz
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
 
end