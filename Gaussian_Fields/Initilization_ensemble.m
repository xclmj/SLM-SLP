%% global generation of permeability ensemble
disp('Reading Gaussian parameter realizations')
tic

nz=1;

load indexMap.mat 

%% properties_clean.gslib contains 10000 laterally correlated random Gaussian 
% properties_clean.gslib contains 101 laterally correlated random Gaussian 
% fields with dimensions 40 x 120 x 20
nx1 = 40; ny1 = 120; nz1 = 20; ng = nx1 * ny1 * nz1; nr = 101;
file = ['properties_clean.gslib'];
fid=fopen(file);
data = fscanf(fid, '%g', nr * ng);
fclose(fid);
data = reshape(data, nr, nx1 * ny1 * nz1);
kmin = min(min(data)); kmax = max(max(data));

data1=[];
for i=1:nr
    data1=[data1 reshape(data(i,:),nx1 * ny1,nz1)];
end

nx =nx1; ny = ny1;
DDx=nx*40; DDy=ny*40; DDz=nz*10;
% set ensemble size ne and create permeability and porosity realizations
np = 300;
jj = 0;
Kperm=[];
Pporo=[];
K1=[];
phi=[];
for ie = 1 : np
    jj = jj + 1;
    kx = data1(:,ie);
    ROCK = 200 * exp(2.5 * kx);
    rock{jj}.perm=log(ROCK);
    rock{jj}.poro= 0.25 * (exp(rock{jj}.perm) ./200).^0.1;
    Kperm(:,jj) = rock{jj}.perm;
    Pporo(:,jj)=rock{jj}.poro;
end

G = cartGrid([nx, ny, nz], [DDx, DDy, DDz]);
G = computeGeometry(G);
Domain_Numx=4; Domain_Numy=5; 

%% Define wells
%Set vertical injectors
iInj = [14,  14,  28, 28, 37, 5, 18];
jInj = [27, 52, 52,  95, 117, 68, 110];

nWI=length(iInj);
q = 200*meter^3/day; bhp = -1;
radius = 0.1;
comp = [1,0];
for k = 1 : numel(iInj)
    name = ['INJ' int2str(k)];
    wells{k} = {iInj(k),jInj(k),q,bhp,radius,comp,name};
end  

% Set vertical producers
iProd = [24, 28,  23, 15,  4, 4 ];
jProd= [14, 35, 68, 90, 12, 118];
nWP=length(iProd);
bhp = 250*barsa(); q = -1;
radius = 0.1;
comp = [0,1];
for k = 1 : numel(iProd)
    name = ['P' int2str(k)];
    wells{nWI+k} = {iProd(k),jProd(k),q,bhp,radius,comp,name};
end
nW=nWI+nWP;         %total number of wells
clear bhp


Initilization_DD;

for ii=1:nz
    
    K=Kperm(1+(ii-1)*nx*ny:ii*nx*ny,1);
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


    for i=1:Domain_Numx-1
        plot([1 ny],[Domain_index(i,4) Domain_index(i,4)],'r--','LineWidth',2);
        hold on
    end

    for i=1:Domain_Numy-1
        plot([Domain_index(1+(i-1)*Domain_Numx,6) Domain_index(1+(i-1)*Domain_Numx,6)],[1 nx],'r--','LineWidth',2);
        hold on
    end
    title(['Permeability'] )

    K=Pporo(1+(ii-1)*nx*ny:ii*nx*ny,1);
    ROCK=nan*ones(nx*ny,1);
    ROCK(indexMap)= Pporo(indexMap,1);
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


    for i=1:Domain_Numx-1
        plot([1 ny],[Domain_index(i,4) Domain_index(i,4)],'r--','LineWidth',2);
        hold on
    end

    for i=1:Domain_Numy-1
        plot([Domain_index(1+(i-1)*Domain_Numx,6) Domain_index(1+(i-1)*Domain_Numx,6)],[1 nx],'r--','LineWidth',2);
        hold on
    end
    title(['Porosity'] )
    
end
