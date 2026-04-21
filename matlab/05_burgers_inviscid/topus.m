function [f] = topus(phi,dx,dt)
N=length(phi);
f=zeros(N,1);
alpha=2;

for j=2:N-1
    ubarra1 = 0.5*(phi(j+1)+phi(j));
    ubarra2 = 0.5*(phi(j)+phi(j-1));
    
    %%%%%%%%%%%%%%%%%%%%%% ubarra1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ubarra1 > 0
        if (phi(j+1)-phi(j-1))==0
            aux1=ubarra1*phi(j);
        else
            cond=(phi(j)-phi(j-1))/(phi(j+1)-phi(j-1));

            if cond>=0 && cond<=1
                aux1=phi(j-1)+(phi(j+1)-phi(j-1))*(alpha*cond^4+(-2*alpha+1)*cond^3+((5*alpha-10)/4)*cond^2+((-alpha+10)/4)*cond);
                aux1 = ubarra1*aux1;
            else
                aux1 = ubarra1*phi(j);
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
                
                if cond>=0 && cond<=1
                    aux1=phi(j+2)+(phi(j)-phi(j+2))*(alpha*cond^4+(-2*alpha+1)*cond^3+((5*alpha-10)/4)*cond^2+((-alpha+10)/4)*cond);
                    aux1 = ubarra1*aux1;
                else
                    aux1 = ubarra1*phi(j+1);
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
                
                if cond>=0 && cond<=1
                    aux2=phi(j-2)+(phi(j)-phi(j-2))*(alpha*cond^4+(-2*alpha+1)*cond^3+((5*alpha-10)/4)*cond^2+((-alpha+10)/4)*cond);
                    aux2 = ubarra2*aux2;
                else
                    aux2 = ubarra2*phi(j-1);
                end
            end          
        end

    else  %ubarra2<=0
        if (phi(j-1)-phi(j+1))==0
            aux2=ubarra2*phi(j);
        else
            cond = (phi(j)-phi(j+1))/(phi(j-1)-phi(j+1));
                
            if cond>=0 && cond<=1
                aux2=phi(j+1)+(phi(j-1)-phi(j+1))*(alpha*cond^4+(-2*alpha+1)*cond^3+((5*alpha-10)/4)*cond^2+((-alpha+10)/4)*cond);
                aux2 = ubarra2*aux2;
            else
                aux2 = ubarra2*phi(j);
            end
        end               
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    conv=0.5*(1/dx)*(aux1-aux2);
    f(j) = phi(j) - (dt)*conv;

    f(1)=phi(1);
    f(N)=phi(N);
end

end