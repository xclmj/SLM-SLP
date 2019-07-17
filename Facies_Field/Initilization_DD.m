clear Domain_Perm 

xx=[];
Domain_sizex=[];
Domain_sizey=[];
yy=0;
if rem(nx,Domain_Numx)==0
   for i=1:Domain_Numx
       Domain_sizex(i)=nx/Domain_Numx;
   end
   
elseif rem(nx,Domain_Numx)>=5
    yy=rem(nx,Domain_Numx);
    for i=1:Domain_Numx
       Domain_sizex(i)=(nx-rem(nx,Domain_Numx))/Domain_Numx;
    end
    Domain_Numx=Domain_Numx+1;
    Domain_sizex(Domain_Numx)=yy;
    
else
    for i=1:Domain_Numx-1
       Domain_sizex(i)=(nx-rem(nx,Domain_Numx))/Domain_Numx;
    end
    Domain_sizex(Domain_Numx)=(nx-rem(nx,Domain_Numx))/Domain_Numx+rem(nx,Domain_Numx);
end


if rem(ny,Domain_Numy)==0
   for i=1:Domain_Numy
       Domain_sizey(i)=ny/Domain_Numy;
   end
   
elseif rem(ny,Domain_Numy)>=5
    yy=rem(ny,Domain_Numy);
    for i=1:Domain_Numy
       Domain_sizey(i)=(ny-rem(ny,Domain_Numy))/Domain_Numy;
    end
    Domain_Numy=Domain_Numy+1;
    Domain_sizey(Domain_Numy)=yy;
    
else
    for i=1:Domain_Numy-1
       Domain_sizey(i)=(ny-rem(ny,Domain_Numy))/Domain_Numy;
    end
    Domain_sizey(Domain_Numy)=(ny-rem(ny,Domain_Numy))/Domain_Numy+rem(ny,Domain_Numy);
end

Domain_index=zeros(Domain_Numx*Domain_Numy,6);
for i=1:Domain_Numy
    for j=1:Domain_Numx
        Domain_index(j+(i-1)*Domain_Numx,1)=j;
        Domain_index(j+(i-1)*Domain_Numx,2)=i;
        
        if j==1
           Domain_index(j+(i-1)*Domain_Numx,3)=1;
           Domain_index(j+(i-1)*Domain_Numx,4)=Domain_sizex(j);
        else
           Domain_index(j+(i-1)*Domain_Numx,3)=sum(Domain_sizex(1:j-1))+1;
           Domain_index(j+(i-1)*Domain_Numx,4)=sum(Domain_sizex(1:j));                
        end
        
        if i==1
           Domain_index(j+(i-1)*Domain_Numx,5)=1;
           Domain_index(j+(i-1)*Domain_Numx,6)=Domain_sizey(i);
        else
           Domain_index(j+(i-1)*Domain_Numx,5)=sum(Domain_sizey(1:i-1))+1;
           Domain_index(j+(i-1)*Domain_Numx,6)=sum(Domain_sizey(1:i));                 
        end
        
    end
end
  
WellI_Domain=zeros(nWI,1);
WellP_Domain=zeros(nWP,1);

for i=1:nWI
    
    for j=1:Domain_Numx*Domain_Numy
            
            if iInj(i)>=Domain_index(j,3) && iInj(i)<=Domain_index(j,4) && jInj(i)>=Domain_index(j,5) && jInj(i)<=Domain_index(j,6)
               
               WellI_Domain(i)=j;
               
            end
            
    end
    
end

for i=1:nWP
    
    for j=1:Domain_Numx*Domain_Numy
            
            if iProd(i)>=Domain_index(j,3) && iProd(i)<=Domain_index(j,4) && jProd(i)>=Domain_index(j,5) && jProd(i)<=Domain_index(j,6)
               
               WellP_Domain(i)=j;
               
            end
            
    end
    
end
