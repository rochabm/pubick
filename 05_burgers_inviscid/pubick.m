function [f] = pubick(phi,dx,dt)
N=length(phi);
f=zeros(N,1);

%---------- parametros para pubick ------------------%
m1=3/10; m2=5/6;
%m1=0.493; m2=0.57;
%-----------------------------------------------------

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
                aux1=phi(j-1)+(phi(j+1)-phi(j-1))*((3/(4*(4*m1-1)^2))*((4*m1)*(4*m1-1)*cond+(2*m1-1)*(2*m1-sqrt(4*m1^2+2*(1-4*m1)*cond))) );
                aux1 = ubarra1*aux1;
            else
                if cond>=0.5 && cond<=1
                    aux1=phi(j-1)+(phi(j+1)-phi(j-1))*(-(1/(4*(4*m2-3)^2))*((4*m2-3)*(5-6*m2+4*(2-3*m2)*cond)+(1-2*m2)*(1-2*m2+sqrt(4*m2^2-2+2*(3-4*m2)*cond))) );
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
                    aux1=phi(j+2)+(phi(j)-phi(j+2))*((3/(4*(4*m1-1)^2))*((4*m1)*(4*m1-1)*cond+(2*m1-1)*(2*m1-sqrt(4*m1^2+2*(1-4*m1)*cond))) );
                    aux1 = ubarra1*aux1;
                else
                    if cond>=0.5 && cond<=1
                        aux1=phi(j+2)+(phi(j)-phi(j+2))*(-(1/(4*(4*m2-3)^2))*((4*m2-3)*(5-6*m2+4*(2-3*m2)*cond)+(1-2*m2)*(1-2*m2+sqrt(4*m2^2-2+2*(3-4*m2)*cond))) );
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
                    aux2=phi(j-2)+(phi(j)-phi(j-2))*((3/(4*(4*m1-1)^2))*((4*m1)*(4*m1-1)*cond+(2*m1-1)*(2*m1-sqrt(4*m1^2+2*(1-4*m1)*cond))) );
                    aux2 = ubarra2*aux2;
                else
                    if cond>=0.5 && cond<=1
                        aux2=phi(j-2)+(phi(j)-phi(j-2))*(-(1/(4*(4*m2-3)^2))*((4*m2-3)*(5-6*m2+4*(2-3*m2)*cond)+(1-2*m2)*(1-2*m2+sqrt(4*m2^2-2+2*(3-4*m2)*cond))) );
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
                aux2=phi(j+1)+(phi(j-1)-phi(j+1))*((3/(4*(4*m1-1)^2))*((4*m1)*(4*m1-1)*cond+(2*m1-1)*(2*m1-sqrt(4*m1^2+2*(1-4*m1)*cond))) );
                aux2 = ubarra2*aux2;
            else
                if cond>=0.5 && cond<=1
                    aux2=phi(j+1)+(phi(j-1)-phi(j+1))*(-(1/(4*(4*m2-3)^2))*((4*m2-3)*(5-6*m2+4*(2-3*m2)*cond)+(1-2*m2)*(1-2*m2+sqrt(4*m2^2-2+2*(3-4*m2)*cond))) );
                    aux2 = ubarra2*aux2;
                else
                    aux2 = ubarra2*phi(j);
                end
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