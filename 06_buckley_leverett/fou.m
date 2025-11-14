function unew = fou(u,dt,dx,tf,c)

k=tf/dt;
N=length(u);
unew=zeros(N,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:round(k)
    for j=2:1:N-1
        if (u(j)==u(j+1))
            v=(dt/dx)*(2*c)*((u(j)-u(j)^2)/(u(j)^2+c*(1-u(j))^2)^2);
        else
            aux1=u(j+1)^2/(u(j+1)^2+c*(1-u(j+1))^2);
            aux2=u(j)^2/(u(j)^2+c*(1-u(j))^2);
            v=(dt/dx)*((aux1-aux2)/(u(j+1)-u(j)));
        end
   
        if (v>0)
            hU=u(j)^2/(u(j)^2+c*(1-u(j))^2);
            top1=hU;          
        else
            hU=u(j+1)^2/(u(j+1)^2+c*(1-u(j+1))^2);
            top1=hU;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (u(j-1)==u(j))
            v=(dt/dx)*(2*c)*((u(j-1)-u(j-1)^2)/(u(j-1)^2+c*(1-u(j-1))^2)^2);
        else
            aux1=(u(j)^2)/(u(j)^2+c*(1-u(j))^2);
            aux2=(u(j-1)^2)/(u(j-1)^2+c*(1-u(j-1))^2);
            v=(dt/dx)*((aux1-aux2)/(u(j)-u(j-1)));
        end
   
        if (v>0)
            hU=(u(j-1)*u(j-1))/(u(j-1)*u(j-1)+c*(1-u(j-1))*(1-u(j-1)));
            top2=hU;
        else
            hU=(u(j)*u(j))/(u(j)*u(j)+c*(1-u(j))*(1-u(j)));
            top2=hU;
        end 
        
        unew(j)=u(j)-(dt/dx)*(top1-top2);  
    end
    unew(1)=u(1);
    unew(N)=u(N);
    u=unew;
end


end