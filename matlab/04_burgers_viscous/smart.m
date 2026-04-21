function [f] = smart(phi,dx,dt,Re)
N=length(phi);
f=zeros(N,1);

for j=2:N-1
    ubarra1 = 0.5*(phi(j+1)+phi(j));
    ubarra2 = 0.5*(phi(j)+phi(j-1));
    
    %%%%%%%%%%%%%%%%%%%%%% ubarra1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ubarra1 > 0
        if (phi(j+1)-phi(j-1))==0
            aux1=ubarra1*phi(j);
        else
            cond=(phi(j)-phi(j-1))/(phi(j+1)-phi(j-1));
            
            if cond>=0 && cond<1/6
                aux1 = ubarra1*(3*phi(j)-2*phi(j-1));             
            else
                if cond>=1/6 && cond<=5/6
                    aux1 = ubarra1*((1/8)*(3*phi(j+1)+6*phi(j)-phi(j-1)));                 
                else
                    if cond>5/6 && cond<=1
                        aux1 = ubarra1*phi(j+1);                  
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
                     
                if cond>=0 && cond<1/6
                    aux1 = ubarra1*(3*phi(j+1)-2*phi(j+2));             
                else
                    if cond>=1/6 && cond<=5/6
                        aux1 = ubarra1*((1/8)*(3*phi(j)+6*phi(j+1)-phi(j+2)));                 
                    else
                        if cond>5/6 && cond<=1
                            aux1 = ubarra1*phi(j);                  
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
            
                if cond>=0 && cond<1/6
                    aux2 = ubarra2*(3*phi(j-1)-2*phi(j-2));             
                else
                    if cond>=1/6 && cond<=5/6
                        aux2 = ubarra2*((1/8)*(3*phi(j)+6*phi(j-1)-phi(j-2)));                 
                    else
                        if cond>5/6 && cond<=1
                            aux2 = ubarra2*phi(j);                  
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
                
            if cond>=0 && cond<1/6
                aux2 = ubarra2*(3*phi(j)-2*phi(j+1));             
            else
                if cond>=1/6 && cond<=5/6
                    aux2 = ubarra2*((1/8)*(3*phi(j-1)+6*phi(j)-phi(j+1)));                 
                else
                    if cond>5/6 && cond<=1
                        aux2 = ubarra2*phi(j-1);                  
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