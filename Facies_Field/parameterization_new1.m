%% parameterization using KL expansion
% for high dimensional model, we use domain decomposition to implement SVD and combine them again.
% Domain_Perm=zeros(max(Domain_sizex)*max(Domain_sizey),np,Domain_Numx*Domain_Numy);
XX=[];
xx=0;
yy=0;
K=[];
clear Layer
for z=1:nz
    
    for k=1:np
        K=Kperm(1+(z-1)*nx*ny:z*nx*ny,k);
        for i=1:Domain_Numx*Domain_Numy   
               XX=[];
               for j=Domain_index(i,5):Domain_index(i,6)
                   XX=[XX; K(Domain_index(i,3)+(j-1)*nx:Domain_index(i,4)+(j-1)*nx)];
               end
               xx=Domain_index(i,1);
               yy=Domain_index(i,2);

               Layer{z}.DomainPerm(1:Domain_sizex(xx)*Domain_sizey(yy),k,i)=XX;
               XX=[];        
        end

    end
    
end

%% step1: generating global basis matrix:
for z=1:nz
    Permeablility=[];         % permeability ensemble
    for i=1:np
        Permeablility=[Permeablility Kperm(1+(z-1)*nx*ny:z*nx*ny,i)]; 
    end

    Layer{z}.mPermeablility=mean(Permeablility,2);
    mmPermeablility=Permeablility-repmat(Layer{z}.mPermeablility,1,np);

    Layer{z}.Cg=mmPermeablility*mmPermeablility'/(np-1);

    [L_G,D_G,R]=svd(Layer{z}.Cg);
    j=sum(D_G);
    psum=0;
    Fala=[];           %%reduced basis for permeability filed
    for i=1:nx*ny
        psum=psum+j(i);
        if psum/sum(j)>EnK_G
           break
        end
    end
    Fala=L_G(:,1:i);
    Layer{z}.fala=Fala*sqrtm(D_G(1:i,1:i));
end

%% step2: generating local basis matrix:
for z=1:nz
    
%    method 1: Domain Decomposition PCA 
    d=0;             %total number of paameter coefficient
    for i=1:Domain_Numx*Domain_Numy
        xx=Domain_index(i,1);
        yy=Domain_index(i,2);
        Permeablility=[];         % permeability ensemble
        Permeablility=Layer{z}.DomainPerm(1:Domain_sizex(xx)*Domain_sizey(yy),:,i);

        Layer{z}.DomainmPerm(1:Domain_sizex(xx)*Domain_sizey(yy),i)=mean(Permeablility,2);
        mmPermeablility=Permeablility-repmat(Layer{z}.DomainmPerm(1:Domain_sizex(xx)*Domain_sizey(yy),i),1,np);

        C=mmPermeablility*mmPermeablility'/(np-1);

        [L1,D1,R]=svd(C);

        Layer{z}.falaK(i)=EnK_num;
        fala11=L1(:,1:EnK_num);
        Layer{z}.Domain_falaK(1:Domain_sizex(xx)*Domain_sizey(yy),1:EnK_num,i)=fala11*sqrtm(D1(1:EnK_num,1:EnK_num));
        d=d+EnK_num;
    end
    Layer{z}.D=d;

end


%% step4: derivative calculation of parameter in global PCA with respect to
%% local PCA
%method1: covariance is non-diagonail
dd3=0;
dd4=0;
TM_GL=[];
for z=1:nz  
    %convert the local basis to global basis
    Basis=zeros(nx*ny,Layer{z}.D);
    for i=1:Domain_Numx*Domain_Numy
        xx=Domain_index(i,1);
        yy=Domain_index(i,2);    
        dd1=sum(Layer{z}.falaK(1:i))-Layer{z}.falaK(i)+1;
        dd2=sum(Layer{z}.falaK(1:i));

        for j=Domain_index(i,5):Domain_index(i,6)
            Basis(Domain_index(i,3)+(j-1)*nx:Domain_index(i,4)+(j-1)*nx,dd1:dd2)=Layer{z}.Domain_falaK(1+(j-Domain_index(i,5))*Domain_sizex(xx):Domain_sizex(xx)+(j-Domain_index(i,5))*Domain_sizex(xx),1:Layer{z}.falaK(i),i);
        end    
    end
    
    Layer{z}.TMGL=(Layer{z}.fala)'*pinv(Layer{z}.Cg)*Basis;
    
    if z==1
        TM_GL(1:size(Layer{z}.fala,2),1:sum(Layer{z}.falaK))= Layer{z}.TMGL  ;
    else
        TM_GL(1+size(Layer{z-1}.fala,2):size(Layer{z-1}.fala,2)+size(Layer{z}.fala,2),1+sum(Layer{z-1}.falaK):sum(Layer{z-1}.falaK)+sum(Layer{z}.falaK))= Layer{z}.TMGL  ;
    end
    
end
