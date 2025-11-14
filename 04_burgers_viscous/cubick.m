function [f] = cubick(phi,dx,dt,Re)
N=length(phi);
f=zeros(N,1);

%------------- parametros para cubick ---------------------%
a=0.25; b=0.45; c=0.256082165462314; d=0.735437227905742;
%a=0.5; b=0.75; c=0.247622138278052; d=0.674286640691196;
%----------------------------------------------------------
A=(3*c-6*a)/(1+3*a-3*c);
B=(3*a)/(1+3*a-3*c);
C=(-1)/(1+3*a-3*c);
p=(3*B-A^2)/3;
q1=(2*A^3-9*A*B)/27;
%----------------------------------------------------------

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
                q=q1+C*cond;
                auxq=sqrt((q^2)/4+(p^3)/27);
                t=(-q/2+auxq)^(1/3)-(q/2+auxq)^(1/3)-A/3;
                auxf=3*b*(t)*(1-t)^2+3*d*(t^2)*(1-t)+t^3;
                aux1=phi(j-1)+(phi(j+1)-phi(j-1))*(auxf);
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
                    q=q1+C*cond;
                    auxq=sqrt((q^2)/4+(p^3)/27);
                    t=(-q/2+auxq)^(1/3)-(q/2+auxq)^(1/3)-A/3;
                    auxf=3*b*(t)*(1-t)^2+3*d*(t^2)*(1-t)+t^3;
                    aux1=phi(j+2)+(phi(j)-phi(j+2))*(auxf);
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
                    q=q1+C*cond;
                    auxq=sqrt((q^2)/4+(p^3)/27);
                    t=(-q/2+auxq)^(1/3)-(q/2+auxq)^(1/3)-A/3;
                    auxf=3*b*(t)*(1-t)^2+3*d*(t^2)*(1-t)+t^3;
                    aux2=phi(j-2)+(phi(j)-phi(j-2))*(auxf);
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
                q=q1+C*cond;
                auxq=sqrt((q^2)/4+(p^3)/27);
                t=(-q/2+auxq)^(1/3)-(q/2+auxq)^(1/3)-A/3;
                auxf=3*b*(t)*(1-t)^2+3*d*(t^2)*(1-t)+t^3;
                aux2=phi(j+1)+(phi(j-1)-phi(j+1))*(auxf);
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