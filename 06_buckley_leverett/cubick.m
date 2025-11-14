function unew=cubick(u,dt,dx,tf,c)
 
k=tf/dt;
N=length(u);
unew=zeros(N,1);
%%%%%%%%%% condiciones pubick %%%%%%%%%%%%%%
%a1=0.25; b1=0.45; c1=0.256082165462314; d1=0.735437227905742;
a1=0.5; b1=0.75; c1=0.247622138278052; d1=0.674286640691196;
A1=(3*c1-6*a1)/(1+3*a1-3*c1);
B1=(3*a1)/(1+3*a1-3*c1); 
C1=(-1)/(1+3*a1-3*c1);
p=(3*B1-A1^2)/3;
q1=(2*A1^3-9*A1*B1)/27;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                        q = q1 + C1*phiU;
                        aux_q = ((q^2)/4 + (p^3)/27)^(1/2); 
                        t = (-q/2 + aux_q)^(1/3) - (q/2 + aux_q)^(1/3) -A1/3; 
                        aux_f = 3*b1*t*(1-t)^2 + 3*d1*(t)^(2)*(1-t) + t^3;
                        top1=hR+(hD-hR)*(aux_f);
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
                        q = q1 + C1*phiU;
                        aux_q = ((q^2)/4 + (p^3)/27)^(1/2); 
                        t = (-q/2 + aux_q)^(1/3) - (q/2 + aux_q)^(1/3) -A1/3; 
                        aux_f = 3*b1*t*(1-t)^2 + 3*d1*(t)^(2)*(1-t) + t^3;
                        top1=hR+(hD-hR)*(aux_f);
                    else
                        top1=hU;
                    end
                end
            else
                top1=hU;
            end
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
            if (j>=3)
                hD=(u(j)*u(j))/(u(j)*u(j)+c*(1-u(j))*(1-u(j)));
                hR=(u(j-2)*u(j-2))/(u(j-2)*u(j-2)+c*(1-u(j-2))*(1-u(j-2)));
                if (hD-hR)==0
                    top2=hU;
                else
                    phiU=(hU-hR)/(hD-hR);
                    if phiU>=0 && phiU<=1
                        q = q1 + C1*phiU;
                        aux_q = ((q^2)/4 + (p^3)/27)^(1/2); 
                        t = (-q/2 + aux_q)^(1/3) - (q/2 + aux_q)^(1/3) -A1/3; 
                        aux_f = 3*b1*t*(1-t)^2 + 3*d1*(t)^(2)*(1-t) + t^3;
                        top2=hR+(hD-hR)*(aux_f);
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
                        q = q1 + C1*phiU;
                        aux_q = ((q^2)/4 + (p^3)/27)^(1/2); 
                        t = (-q/2 + aux_q)^(1/3) - (q/2 + aux_q)^(1/3) -A1/3; 
                        aux_f = 3*b1*t*(1-t)^2 + 3*d1*(t)^(2)*(1-t) + t^3;
                        top2=hR+(hD-hR)*(aux_f);
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





