function unew=fdhpus(u,dt,dx,tf,c)
 
k=tf/dt;
N=length(u);
unew=zeros(N,1);
%%%%%%%%%%%%%%%%%%%
a=1.5;
b=0;
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
            if (j>=2)
                hD=u(j+1)^2/(u(j+1)^2+c*(1-u(j+1))^2);
                hR=u(j-1)^2/(u(j-1)^2+c*(1-u(j-1))^2);
                if ((hD-hR)==0)
                    top1=hU;
                else
                    phiU=(hU-hR)/(hD-hR);
                    if phiU>=0 && phiU<=1
                        top1=hR+(hD-hR)*(4*(a+b-3)*phiU^5-2*(6*a+4*b-17)*phiU^4+(13*a+5*b-34)*phiU^3-(6*a+b-13)*phiU^2+a*phiU);
                    else
                        top1=hU;
                    end
                end
            else
                top1=hU;
            end          
        else
            hU=u(j+1)^2/(u(j+1)^2+c*(1-u(j+1))^2);
            if (j<=N-2)
                hD=(u(j)^2)/(u(j)^2+c*(1-u(j))^2);
                hR=(u(j+2)^2)/(u(j+2)^2+c*(1-u(j+2))^2);

                if ((hD-hR)==0)
                    top1=hU;
                else        
                    phiU=(hU-hR)/(hD-hR);
                    if phiU>=0 && phiU<=1
                        top1=hR+(hD-hR)*(4*(a+b-3)*phiU^5-2*(6*a+4*b-17)*phiU^4+(13*a+5*b-34)*phiU^3-(6*a+b-13)*phiU^2+a*phiU);
                    else
                        top1=hU;
                    end
                end
            else
                top1=hU;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (u(j-1)==u(j))
            v=(dt/dx)*(2*c)*((u(j-1)-u(j-1)^2)/(u(j-1)^2+c*(1-u(j-1))^2)^2);
        else
            aux1=(u(j)^2)/(u(j)^2+c*(1-u(j))^2);
            aux2=(u(j-1)^2)/(u(j-1)^2+c*(1-u(j-1))^2);
            v=(dt/dx)*((aux1-aux2)/(u(j)-u(j-1)));
        end
   
        if (v>0)
            hU=(u(j-1)*u(j-1))/(u(j-1)*u(j-1)+c*(1-u(j-1))*(1-u(j-1)));
            if (j>=3)
                hD=(u(j)*u(j))/(u(j)*u(j)+c*(1-u(j))*(1-u(j)));
                hR=(u(j-2)*u(j-2))/(u(j-2)*u(j-2)+c*(1-u(j-2))*(1-u(j-2)));
                if (hD-hR)==0
                    top2=hU;
                else
                    phiU=(hU-hR)/(hD-hR);
                    if phiU>=0 && phiU<=1
                        top2=hR+(hD-hR)*(4*(a+b-3)*phiU^5-2*(6*a+4*b-17)*phiU^4+(13*a+5*b-34)*phiU^3-(6*a+b-13)*phiU^2+a*phiU);
                    else
                        top2=hU;
                    end
                end
            else
                top2=hU;
            end
        else
            hU=(u(j)*u(j))/(u(j)*u(j)+c*(1-u(j))*(1-u(j)));
            if (j<=N-2)
                hD=(u(j-1)*u(j-1))/(u(j-1)*u(j-1)+c*(1-u(j-1))*(1-u(j-1)));
                hR=(u(j+1)*u(j+1))/(u(j+1)*u(j+1)+c*(1-u(j+1))*(1-u(j+1)));
                if (hD-hR)==0
                    top2=hU;
                else
                    phiU=(hU-hR)/(hD-hR);
                    if phiU>=0 && phiU<=1
                        top2=hR+(hD-hR)*(4*(a+b-3)*phiU^5-2*(6*a+4*b-17)*phiU^4+(13*a+5*b-34)*phiU^3-(6*a+b-13)*phiU^2+a*phiU);
                    else
                        top2=hU;
                    end
                end
            else
                top2=hU;
            end
        end 
        
        unew(j)=u(j)-(dt/dx)*(top1-top2);  
    end
    unew(1)=u(1);
    unew(N)=u(N);
    u=unew;
end



end





