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

EnK_G=0.95;
EnK_num=5;
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
    Layer{z}.L=L_G*sqrtm(D_G);
end

variance=j;
Energy = cumsum(variance) / sum(variance); % relative energy
%% Plot variance and relative energy
figure;
subplot(1,2,1)
hold on
semilogy(1:length(variance),variance, 'b','linewidth',2);
xlabel('Number of retained eigenpairs','FontName','Times New Roman','FontSize',14,'FontWeight','n')
ylabel('Variance (eigenvalue)','FontName','Times New Roman','FontSize',14,'FontWeight','n')

%  Relative energy
subplot(1,2,2)
hold on
semilogx(1:length(Energy),Energy, 'b','linewidth',2);
xlabel('Number of retained eigenpairs','FontName','Times New Roman','FontSize',14,'FontWeight','n')
ylabel('Relative energy','FontName','Times New Roman','FontSize',14,'FontWeight','n')

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

%% step3: optimization-based local parameterization technique
%method1: covariance is non-diagonail
K_LPCA=[];
K_OLPCA=[];
K_GPCA=[];
for z=1:nz
    
    n=0;
    for jj=1:1
        verify=rand(Layer{z}.D,1)-0.5*ones(Layer{z}.D,1);
        verify=2*verify;
        KK=[];
        for i=1:Domain_Numx*Domain_Numy
            xx=Domain_index(i,1);
            yy=Domain_index(i,2);    
            dd1=sum(Layer{z}.falaK(1:i))-Layer{z}.falaK(i)+1;
            dd2=sum(Layer{z}.falaK(1:i));
            Layer{z}.DomainPerm(1:Domain_sizex(xx)*Domain_sizey(yy),1,i)=Layer{z}.DomainmPerm(1:Domain_sizex(xx)*Domain_sizey(yy),i)+Layer{z}.Domain_falaK(1:Domain_sizex(xx)*Domain_sizey(yy),1:Layer{z}.falaK(i),i)*verify(dd1:dd2);

            k=(z-1)*nx*ny;
            for j=Domain_index(i,5):Domain_index(i,6)
                KK(Domain_index(i,3)+(j-1)*nx+k:Domain_index(i,4)+(j-1)*nx+k,1)=Layer{z}.DomainPerm(1+(j-Domain_index(i,5))*Domain_sizex(xx):Domain_sizex(xx)+(j-Domain_index(i,5))*Domain_sizex(xx),1,i);
            end    
        end 
    end

    K=KK(1+(z-1)*nx*ny:z*nx*ny);
    K_LPCA=[K_LPCA ; K];    
    K=K-Layer{z}.mPermeablility;
   
    
    verify_LPCA= (Layer{z}.fala)'*pinv(Layer{z}.Cg)*K;
    K_OLPCA=[K_OLPCA ; Layer{z}.mPermeablility+Layer{z}.fala*verify_LPCA];

    verify_G=(rand(nx*ny,1)-0.5*ones(nx*ny,1))*2;
    verify_G(1:length(verify_LPCA))=verify_LPCA;
    K_GPCA=[K_GPCA ; Layer{z}.mPermeablility+Layer{z}.L*verify_G];
    

end

%% step4: derivative calculation of parameter in global PCA with respect to
%% local PCA
%method1: covariance is non-diagonail
dd3=0;
dd4=0;
TM_GL=[];
x=0;
y=0;
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
        x=x+size(Layer{z}.fala,2);
        y=y+sum(Layer{z}.falaK);
        TM_GL(1:size(Layer{z}.fala,2),1:sum(Layer{z}.falaK))= Layer{z}.TMGL  ;
    else
        x=x+size(Layer{z}.fala,2);
        y=y+sum(Layer{z}.falaK);
        TM_GL(x-size(Layer{z}.fala,2)+1:x,y-sum(Layer{z}.falaK)+1:y)= Layer{z}.TMGL  ;
    end

end

for i=1:nz
    
    K=K_LPCA(1+(i-1)*nx*ny:i*nx*ny,1);
    ROCK=nan*ones(nx*ny,1);
    ROCK(indexMap)= K(indexMap);
    ROCK=reshape(ROCK,nx,ny,1);
    figure('color',[1,1,1])
    cm = colormap(jet);
    gca=pcolor(ROCK); colorbar; box on
    set(gca, 'LineStyle','none');
    hold on
    plot(jInj, iInj, '^k', 'MarkerFaceColor', [0.99 0.99 0.99], 'MarkerSize', 8); % used to differentiate injection and production wells
    for index=1:length(jInj)
        text(jInj(index)-1, iInj(index)-2, ['I' num2str(index)], 'FontName', 'Times New Roman', 'color', 'k', 'FontSize', 12, 'FontWeight', 'b');
    end

    plot(jProd, iProd, 'ok', 'MarkerFaceColor',[0.99 0.99 0.99], 'MarkerSize',8);
    for index=1:length(jProd)
        text(jProd(index)-2, iProd(index)-2, ['P' num2str(index)], 'FontName', 'Times New Roman', 'color', 'k', 'FontSize', 12, 'FontWeight', 'b'); 
    end
    title(['LP'])

    K=K_OLPCA(1+(i-1)*nx*ny:i*nx*ny,1);
    ROCK=nan*ones(nx*ny,1);
    ROCK(indexMap)= K(indexMap);
    ROCK=reshape(ROCK,nx,ny,1);
    figure('color',[1,1,1])
    cm = colormap(jet);
    gca=pcolor(ROCK); colorbar; box on
    set(gca, 'LineStyle','none');
    hold on
    plot(jInj, iInj, '^k', 'MarkerFaceColor', [0.99 0.99 0.99], 'MarkerSize', 8); % used to differentiate injection and production wells
    for index=1:length(jInj)
        text(jInj(index)-1, iInj(index)-2, ['I' num2str(index)], 'FontName', 'Times New Roman', 'color', 'k', 'FontSize', 12, 'FontWeight', 'b');
    end

    plot(jProd, iProd, 'ok', 'MarkerFaceColor',[0.99 0.99 0.99], 'MarkerSize',8);
    for index=1:length(jProd)
        text(jProd(index)-2, iProd(index)-2, ['P' num2str(index)], 'FontName', 'Times New Roman', 'color', 'k', 'FontSize', 12, 'FontWeight', 'b'); 
    end
    title(['SLP'])

end

