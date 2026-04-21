function unew=adbquickest(u,dt,dx,tf,c,CFL)
 
k=tf/dt;
N=length(u);
unew=zeros(N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
numa = ((1/3)-(1/2)*abs(CFL)+(1/6)*CFL^2);
dena = (2*((7/12)-(1/2)*CFL-(1/4)*abs(CFL)+(1/6)*CFL^2));
a    = numa/dena;

numb = ((1/6)*CFL^2-(2/3)-(1/2)*abs(CFL)+CFL);
denb = (2*((-5/12)-(1/4)*abs(CFL)+(1/6)*CFL^2+(1/2)*CFL));
b    = numb/denb;
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
            
                hD=u(j+1)^2/(u(j+1)^2+c*(1-u(j+1))^2);
                hR=u(j-1)^2/(u(j-1)^2+c*(1-u(j-1))^2);

                if ((hD-hR)==0)
                    top1=hU;
                else
                    phiU=(hU-hR)/(hD-hR);
                    A=hU-hR;
                    B=hD-hR;
                    C=hR;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if phiU>=a && phiU<=b
                        top1 =C+A+0.5*(1-abs(CFL))*(B-A)-(1/6)*(1-CFL^2)*(B-2*A);
                    else
                        if phiU>=0 && phiU<a
                            top1=C+(2-CFL)*A;
                        else
                            if phiU>b && phiU<=1
                                top1=C+B-B*CFL+A*CFL;
                            else
                                top1=hU;
                            end
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                    A=hU-hR;
                    B=hD-hR;
                    C=hR;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if phiU>=a && phiU<=b
                        top1 =C+A+0.5*(1-abs(CFL))*(B-A)-(1/6)*(1-CFL^2)*(B-2*A);
                    else
                        if phiU>=0 && phiU<a
                            top1=C+(2-CFL)*A;
                        else
                            if phiU>b && phiU<=1
                                top1=C+B-B*CFL+A*CFL;
                            else
                                top1=hU;
                            end
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                    A=hU-hR;
                    B=hD-hR;
                    C=hR;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if phiU>=a && phiU<=b
                        top2 =C+A+0.5*(1-abs(CFL))*(B-A)-(1/6)*(1-CFL^2)*(B-2*A);
                    else
                        if phiU>=0 && phiU<a
                            top2=C+(2-CFL)*A;
                        else
                            if phiU>b && phiU<=1
                                top2=C+B-B*CFL+A*CFL;
                            else
                                top2=hU;
                            end
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            else
                top2=hU;
            end
        else
            hU=(u(j)*u(j))/(u(j)*u(j)+c*(1-u(j))*(1-u(j)));
            
                hD=(u(j-1)*u(j-1))/(u(j-1)*u(j-1)+c*(1-u(j-1))*(1-u(j-1)));
                hR=(u(j+1)*u(j+1))/(u(j+1)*u(j+1)+c*(1-u(j+1))*(1-u(j+1)));
                if (hD-hR)==0
                    top2=hU;
                else
                    phiU=(hU-hR)/(hD-hR);              
                    A=hU-hR;
                    B=hD-hR;
                    C=hR;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if phiU>=a && phiU<=b
                        top2 =C+A+0.5*(1-abs(CFL))*(B-A)-(1/6)*(1-CFL^2)*(B-2*A);
                    else
                        if phiU>=0 && phiU<a
                            top2=C+(2-CFL)*A;
                        else
                            if phiU>b && phiU<=1
                                top2=C+B-B*CFL+A*CFL;
                            else
                                top2=hU;
                            end
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end

        end 
       
        
        unew(j)=u(j)-(dt/dx)*(top1-top2);  
    end
    unew(1)=u(1);
    unew(N)=u(N);
    u=unew;
end



end
