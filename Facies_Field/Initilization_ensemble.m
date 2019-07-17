%% parameterization facies using KL expansion
disp('Reading facies realizations')
clear wells
ROCK=[];

X = load('S22_1000CondRealsTo16Wells_60x60_3mudHardData.txt'); % model: 60x60, binary, 2D, conditioned to hard data, 16 well
[Nc, np] = size(X);% Nc: number of cells; Nr: number of realizations

np=100;
nx=40; ny=40; nz=1;
K_sand=2000;  K_mud=20;
C=reshape(X(:,1:np),60,60,np);
C=C(1:nx,1:ny,1:np);
X=reshape(C,nx*ny,np);
Kperm=X;

D=zeros(nx*ny,np);
for i=1:np
    for j=1:nx*ny
        if X(j,i)==1;
           D(j,i)=2000;
        else
           D(j,i)=20;
        end
    end
    rock{i}.perm=X(:,i);
    rock{i}.poro = 0.25 * (D(:,i) ./200).^0.1;
end    

DDx=nx*40; DDy=ny*40; DDz=10;
G = cartGrid([nx, ny, nz], [DDx, DDy, DDz]);
G = computeGeometry(G);    

iInj = [5, 5, 33, 33]; % 4 wells configuration, injection wells
jInj = [5, 33, 6, 37]; % 4 wells configuration, injection wells
q = 200*meter^3/day; bhp = -1;
radius = 0.1;
comp = [1,0];
for k = 1 : numel(iInj)
    name = ['I' int2str(k)];
    wells{k} = {iInj(k),jInj(k),q,bhp,radius,comp,name};
end
nWI = length(iInj);

iProd = [13, 13, 20, 13, 25, 28 ]; % 16 wells configuration, production wells
jProd = [ 7, 21, 35, 31, 7, 22]; % 16 wells configuration, production wells
for i=1:length(iProd)
    bhp(i) = (250-i*0)*barsa(); q(i) = -1;
end

radius = 0.1;
comp = [0,1];
for k = 1 : numel(iProd)
    name = ['P' int2str(k)];
    wells{nWI+k} = {iProd(k),jProd(k),q(k),bhp(k),radius,comp,name};
end
nWP = length(iProd);
nW=nWI+nWP;
clear bhp 

Domain_Numx=3; Domain_Numy=3;
    
Initilization_DD;

for i=1:nz
    
    K=rock{randi(np,1)}.perm(1+(i-1)*nx*ny:i*nx*ny,1);
    ROCK=nan*ones(nx*ny,1);
    ROCK= K;
    ROCK=reshape(ROCK,nx,ny,1);
    figure('color',[1,1,1])
    cm = colormap(jet);
    gca=pcolor(ROCK); 
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
    title([num2str(i) 'st layer'])

    for i=1:Domain_Numy-1
        plot([1 nx], [Domain_index(1+(i-1)*Domain_Numx,6) Domain_index(1+(i-1)*Domain_Numx,6)],'r--','LineWidth',2);
        hold on
    end

    for i=1:Domain_Numx-1
        plot([Domain_index(i,4) Domain_index(i,4)], [1 ny],'r--','LineWidth',2);
        hold on
    end

end

