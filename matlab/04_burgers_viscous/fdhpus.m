function [f] = fdhpus(phi,dx,dt,Re)
N=length(phi);
f=zeros(N,1);

%--------- parametros fdhpus ------------
a=1.5; b=0; % a=\theta_1, b=\theta_2 
%a=2; b=0; %(wei)
%----------------------------------------

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
                aux1=phi(j-1)+(phi(j+1)-phi(j-1))*(4*(a+b-3)*cond^5-2*(6*a+4*b-17)*cond^4+(13*a+5*b-34)*cond^3-(6*a+b-13)*cond^2+a*cond);
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
                    aux1=phi(j+2)+(phi(j)-phi(j+2))*(4*(a+b-3)*cond^5-2*(6*a+4*b-17)*cond^4+(13*a+5*b-34)*cond^3-(6*a+b-13)*cond^2+a*cond);
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
                    aux2=phi(j-2)+(phi(j)-phi(j-2))*(4*(a+b-3)*cond^5-2*(6*a+4*b-17)*cond^4+(13*a+5*b-34)*cond^3-(6*a+b-13)*cond^2+a*cond);
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
                aux2=phi(j+1)+(phi(j-1)-phi(j+1))*(4*(a+b-3)*cond^5-2*(6*a+4*b-17)*cond^4+(13*a+5*b-34)*cond^3-(6*a+b-13)*cond^2+a*cond);
                aux2 = ubarra2*aux2;
            else
                aux2 = ubarra2*phi(j);
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