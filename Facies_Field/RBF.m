function [BB]=RBF(j,falaK,Domain_index,Domain_Numx,Domain_Numy,nz)
 
% identify the neighboring subdomain
if Domain_index(j,1)==1 && Domain_index(j,2)==1
   BB=[];

   for z=1:nz
       dd1=sum(falaK(1:j))-falaK(j)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:j))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j+1;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j+Domain_Numx;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];   

       k=j+Domain_Numx+1;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];  
    end



elseif Domain_index(j,1)==1 && Domain_index(j,2)==Domain_Numy

   BB=[];

   for z=1:nz
       k=j-Domain_Numx;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j-Domain_Numx+1;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];   

       dd1=sum(falaK(1:j))-falaK(j)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:j))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j+1;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];
    end


elseif Domain_index(j,1)==Domain_Numx && Domain_index(j,2)==1

   BB=[];
   for z=1:nz 
       k=j-1;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];   

       k=j+Domain_Numx-1;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j+Domain_Numx;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2]; 
  end

elseif Domain_index(j,1)==Domain_Numx && Domain_index(j,2)==Domain_Numy

   BB=[];
   for z=1:nz
       k=j-Domain_Numx-1;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j-Domain_Numx;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];   

       k=j-1;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];  
    end

 elseif Domain_index(j,1)==1 && Domain_index(j,2)>1

   BB=[];
   for z=1:nz
       k=j-Domain_Numx;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];   

       k=j+1;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j+Domain_Numx;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];
   end

 elseif Domain_index(j,1)==Domain_Numx && Domain_index(j,2)>1

   BB=[];
   for z=1:nz
       k=j-Domain_Numx;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j-1;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];   

       k=j;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j+Domain_Numx;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];
    end

 elseif Domain_index(j,1)>1 && Domain_index(j,2)==1

   BB=[];
   for z=1:nz
       k=j-1;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];   

       k=j+1;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j+Domain_Numx;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];
    end

 elseif Domain_index(j,1)>1 && Domain_index(j,2)==Domain_Numy

   BB=[];
   for z=1:nz
       k=j-Domain_Numx;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j-1;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];   

       k=j;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j+1;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];
    end

else

   BB=[];
   for z=1:nz
       k=j-Domain_Numx;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j-1;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];   

       k=j;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j+1;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];

       k=j+Domain_Numx;
       dd1=sum(falaK(1:k))-falaK(k)+1+(z-1)*sum(falaK);
       dd2=sum(falaK(1:k))+(z-1)*sum(falaK);
       BB=[BB dd1:dd2];   
    end

end
        
        

