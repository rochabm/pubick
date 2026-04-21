function [f] = pubick_calc2(phi,dx,dt,Re)
% Esquema PUBICK  con los caculos realizados para los valores %m1=3/10 y m2=5/6;
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

            if cond>=0 && cond<0.5
                aux1=phi(j-1)+(phi(j+1)-phi(j-1))*(1.521604938271605*cond + 0.011113651374282*(0.972196 - 1.944*cond)^(1/2)- 0.010958060255042);
                aux1 = ubarra1*aux1;
            else
                if cond>=0.5 && cond<=1
                    aux1=phi(j-1)+(phi(j+1)-phi(j-1))*((29*cond)/72 + (175*((36*cond)/25 - 1751/2500)^(1/2))/2592 + 2795/5184);
                    aux1 = ubarra1*aux1;
                else
                    aux1 = ubarra1*phi(j);
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
                
                if cond>=0 && cond<0.5
                    aux1=phi(j+2)+(phi(j)-phi(j+2))*(1.521604938271605*cond + 0.011113651374282*(0.972196 - 1.944*cond)^(1/2)- 0.010958060255042);
                    aux1 = ubarra1*aux1;
                else
                    if cond>=0.5 && cond<=1
                        aux1=phi(j+2)+(phi(j)-phi(j+2))*((29*cond)/72 + (175*((36*cond)/25 - 1751/2500)^(1/2))/2592 + 2795/5184);
                        aux1 = ubarra1*aux1;
                    else
                        aux1 = ubarra1*phi(j+1);
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
                
                if cond>=0 && cond<0.5
                    aux2=phi(j-2)+(phi(j)-phi(j-2))*(1.521604938271605*cond + 0.011113651374282*(0.972196 - 1.944*cond)^(1/2)- 0.010958060255042);
                    aux2 = ubarra2*aux2;
                else
                    if cond>=0.5 && cond<=1
                        aux2=phi(j-2)+(phi(j)-phi(j-2))*((29*cond)/72 + (175*((36*cond)/25 - 1751/2500)^(1/2))/2592 + 2795/5184);
                        aux2 = ubarra2*aux2;
                    else
                        aux2 = ubarra2*phi(j-1);
                    end
                end
            end          
        end

    else  %ubarra2<=0
        if (phi(j-1)-phi(j+1))==0
            aux2=ubarra2*phi(j);
        else
            cond = (phi(j)-phi(j+1))/(phi(j-1)-phi(j+1));
            
            if cond>=0 && cond<0.5
                aux2=phi(j+1)+(phi(j-1)-phi(j+1))*(1.521604938271605*cond + 0.011113651374282*(0.972196 - 1.944*cond)^(1/2)- 0.010958060255042);
                aux2 = ubarra2*aux2;
            else
                if cond>=0.5 && cond<=1
                    aux2=phi(j+1)+(phi(j-1)-phi(j+1))*((29*cond)/72 + (175*((36*cond)/25 - 1751/2500)^(1/2))/2592 + 2795/5184);
                    aux2 = ubarra2*aux2;
                else
                    aux2 = ubarra2*phi(j);
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