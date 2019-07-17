figure(1)
subplot(1,2,1)
boxplot(RMSE_ensemble_EnRML')
subplot(1,2,2)
boxplot(RMSE_K_EnRML')

figure(2)
subplot(1,2,1)
boxplot(RMSE_ensemble_ESMDA')
subplot(1,2,2)
boxplot(RMSE_K_ESMDA')

rock{1}.perm=exp(K_ref);
rock{1}.poro = 0.25 * (rock{1}.perm ./200).^0.1;
rock{1}.perm=convertFrom(rock{1}.perm,milli*darcy);

Num_RML=10;
for jj=1:Num_RML

    verify2=verify_RML_Init(:,jj);

    for z=1:nz
        xx=0;
        for zz=1:z
            xx=xx+Layer{zz}.D;
        end
        verify1=verify2(xx-Layer{zz}.D+1:xx);
        K=Layer{z}.mPermeablility+Layer{z}.fala*Layer{z}.TMGL*verify1;
        PORO=0.25 * (exp(K) ./200).^0.1;
        ROCK=convertFrom(exp(K),milli*darcy);
        rock{jj+1}.perm(1+(z-1)*nx*ny:z*nx*ny,1)=ROCK;
        rock{jj+1}.poro(1+(z-1)*nx*ny:z*nx*ny,1)=PORO;   

    end

end

for jj=1:Num_RML

    verify2=verify_ESMDA(:,jj);

    for z=1:nz
        xx=0;
        for zz=1:z
            xx=xx+Layer{zz}.D;
        end
        verify1=verify2(xx-Layer{zz}.D+1:xx);
        K=Layer{z}.mPermeablility+Layer{z}.fala*Layer{z}.TMGL*verify1;
        PORO=0.25 * (exp(K) ./200).^0.1;
        ROCK=convertFrom(exp(K),milli*darcy);
        rock{jj+1+Num_RML}.perm(1+(z-1)*nx*ny:z*nx*ny,1)=ROCK;
        rock{jj+1+Num_RML}.poro(1+(z-1)*nx*ny:z*nx*ny,1)=PORO;   

    end

end

tic
Prediction_EDS
toc

TIME0=1:1:istep;
TIME0=TIME0*dT/3600/24/365;

figure(10)
for j=1:round(nWP)
    subplot(1,round(nWP),j);
    plot(TIME0,wellDataSGSPC(nWI+j).flx(1,:),'b','LineWidth',2);
    hold on
    for i=1:Num_RML
        plot(TIME0,wellDataSGSPC(nWI+j).flx(1+i,:),'color',[0.5 0.5 0.5]);
    end
    for i=1:Num_RML
        plot(TIME0,wellDataSGSPC(nWI+j).flx(1+Num_RML+i,:),'r');
    end    
end
xlabel('production time, year')
ylabel('Fluid Rate, m3/d')


figure(12)
for j=1:round(nWP)
    subplot(1,round(nWP),j);
    plot(TIME0,wellDataSGSPC(nWI+j).ffl(1,:),'b','LineWidth',2);
    hold on
    for i=1:Num_RML
        plot(TIME0,wellDataSGSPC(nWI+j).ffl(1+i,:),'color',[0.5 0.5 0.5]);
    end
    for i=1:Num_RML
        plot(TIME0,wellDataSGSPC(nWI+j).ffl(1+Num_RML+i,:),'r');
    end   
end 
xlabel('production time, year')
ylabel('Wtercut, m3/d')


ROCK=nan*ones(nx*ny,1);
ROCK= K_ref;
ROCK=reshape(ROCK,nx,ny,1);
figure('color',[1,1,1])
cm = colormap(jet);
gca=pcolor(ROCK); colorbar; box on;caxis([1, 8])
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

for jj=1:5

    verify2=verify_RML_Init(:,jj);

    for z=1:nz
        xx=0;
        for zz=1:z
            xx=xx+Layer{zz}.D;
        end
        verify1=verify2(xx-Layer{zz}.D+1:xx);
        K=Layer{z}.mPermeablility+Layer{z}.fala*Layer{z}.TMGL*verify1;
    end

    ROCK=nan*ones(nx*ny,1);
    ROCK= K;
    ROCK=reshape(ROCK,nx,ny,1);
    figure('color',[1,1,1])
    cm = colormap(jet);
    gca=pcolor(ROCK); colorbar; box on;caxis([1, 8])
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


end
   
for jj=1:5

    verify2=verify_ESMDA(:,jj);

    for z=1:nz
        xx=0;
        for zz=1:z
            xx=xx+Layer{zz}.D;
        end
        verify1=verify2(xx-Layer{zz}.D+1:xx);
        K=Layer{z}.mPermeablility+Layer{z}.fala*Layer{z}.TMGL*verify1;
    end

    ROCK=nan*ones(nx*ny,1);
    ROCK= K;
    ROCK=reshape(ROCK,nx,ny,1);
    figure('color',[1,1,1])
    cm = colormap(jet);
    gca=pcolor(ROCK); colorbar; box on;caxis([1, 8])
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


end

