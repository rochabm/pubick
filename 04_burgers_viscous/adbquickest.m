function [f] = adbquickest(phi,dx,dt,Re,cfl)
N=length(phi);
f=zeros(N,1);

%---------------------- parametros -----------------------
alfad  = (1/2)*(1-abs(cfl))-(1/6)*(1-cfl^2);
alfau  = 1-(1/2)*(1-abs(cfl))+(1/3)*(1-cfl^2);
alfar  = (1/6)*(1-cfl^2);

numa = (2-(3*abs(cfl))+(cfl^2));
dena = (7-(6*cfl)-(3*abs(cfl))+(2*cfl^2));
a    = numa/dena;

numb = ((cfl^2)-4-(3*abs(cfl))+(6*cfl));
denb = (-5-(3*abs(cfl))+(2*cfl^2)+(6*cfl));
b    = numb/denb;
%---------------------------------------------------------

for j=2:N-1
    ubarra1 = 0.5*(phi(j+1)+phi(j));
    ubarra2 = 0.5*(phi(j)+phi(j-1));
    
    %%%%%%%%%%%%%%%%%%%%%% ubarra1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ubarra1 > 0
        if (phi(j+1)-phi(j-1))==0
            aux1=ubarra1*phi(j);
        else
            cond=(phi(j)-phi(j-1))/(phi(j+1)-phi(j-1));
               
            if 0<=cond && cond<a
                aux1 = ubarra1*((2-cfl)*phi(j)-(1-cfl)*phi(j-1));
            else
                if a<=cond && cond<=b
                    aux1 = ubarra1*(alfad*phi(j+1)+alfau*phi(j)-alfar*phi(j-1));
                else
                    if b<cond && cond<=1
                        aux1 = ubarra1*((1.0-cfl)*phi(j+1)+cfl*phi(j));
                    else
                        aux1 = ubarra1*phi(j);
                    end
                end
            end
        end

    else  %ubarra1<=0
        if j==N-1
            aux1 = ubarra1*0.5*(phi(j)+phi(j+1));
        else
            if (phi(j)-phi(j+2))==0
                aux1=ubarra1*phi(j+1);
            else
                cond = (phi(j+1)-phi(j+2))/(phi(j)-phi(j+2));
                
                if 0<=cond && cond<a
                    aux1 = ubarra1*((2-cfl)*phi(j+1)-(1-cfl)*phi(j+2));
                else
                    if a<=cond && cond<=b
                        aux1 = ubarra1*(alfad*phi(j)+alfau*phi(j+1)-alfar*phi(j+2));
                    else
                        if b<cond && cond<=1
                            aux1 = ubarra1*((1.0-cfl)*phi(j)+cfl*phi(j+1));
                        else
                            aux1 = ubarra1*phi(j+1);
                        end
                    end
                end
            end          
        end
    end

    %%%%%%%%%%%%%%%%%%%%%% ubarra2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ubarra2 > 0
        if j==2
            aux2 = ubarra2*0.5*(phi(j)+phi(j-1));
        else
            if (phi(j)-phi(j-2))==0
                aux2=ubarra2*phi(j-1);
            else
                cond = (phi(j-1)-phi(j-2))/(phi(j)-phi(j-2));
                
                if 0<=cond && cond<a
                    aux2 = ubarra2*((2-cfl)*phi(j-1)-(1-cfl)*phi(j-2));
                else
                    if a<=cond && cond<=b
                        aux2 = ubarra2*(alfad*phi(j)+alfau*phi(j-1)-alfar*phi(j-2));
                    else
                        if b<cond && cond<=1
                            aux2 = ubarra2*((1.0-cfl)*phi(j)+cfl*phi(j-1));
                        else
                            aux2 = ubarra2*phi(j-1);
                        end
                    end
                end
            end          
        end

    else  %ubarra2<=0
        if (phi(j-1)-phi(j+1))==0
            aux2=ubarra2*phi(j);
        else
            cond = (phi(j)-phi(j+1))/(phi(j-1)-phi(j+1));

             if 0<=cond && cond<a
                 aux2 = ubarra2*((2-cfl)*phi(j)-(1-cfl)*phi(j+1));
             else
                 if a<=cond && cond<=b
                     aux2 = ubarra2*(alfad*phi(j-1)+alfau*phi(j)-alfar*phi(j+1));
                 else
                     if b<cond && cond<=1
                         aux2 = ubarra2*((1.0-cfl)*phi(j-1)+cfl*phi(j));
                     else
                         aux2 = ubarra2*phi(j);
                     end
                 end
             end
        end               
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    conv=0.5*(1/dx)*(aux1-aux2);
    f(j) = phi(j) - (dt)*conv + ((dt)/(Re*dx^2))*(phi(j+1)-2*phi(j)+phi(j-1));

    f(1)=phi(1);
    f(N)=phi(N);
end

end