function [f] = pubick(phi)
N=length(phi);
f=zeros(N,1);
m1=3/10; m2=5/6;
%m1=0.493; m2=0.57;

for j=2:N-1
    if (phi(j+1)-phi(j-1))==0
        f(j)=phi(j);
    else    
        cond=(phi(j)-phi(j-1))/(phi(j+1)-phi(j-1));
            
        if cond>=0 && cond<0.5
            f(j)=phi(j-1)+(phi(j+1)-phi(j-1))*((3/(4*(4*m1-1)^2))*((4*m1)*(4*m1-1)*cond+(2*m1-1)*(2*m1-sqrt(4*m1^2+2*(1-4*m1)*cond))) );
  
        else
            if cond>=0.5 && cond<=1
                f(j)=phi(j-1)+(phi(j+1)-phi(j-1))*(-(1/(4*(4*m2-3)^2))*((4*m2-3)*(5-6*m2+4*(2-3*m2)*cond)+(1-2*m2)*(1-2*m2+sqrt(4*m2^2-2+2*(3-4*m2)*cond))) );
            else
                f(j)=phi(j);
            end
        end
    end
end   
f(1)=f(N-1);
f(N)=f(2);
end

