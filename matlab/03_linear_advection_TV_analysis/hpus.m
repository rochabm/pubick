function [f] = hpus(phi)

N=length(phi);
f=zeros(N);

for j=2:N-1
    if (phi(j+1)-phi(j-1))==0
        f(j)=phi(j);
    else
        cond=(phi(j)-phi(j-1))/(phi(j+1)-phi(j-1));
            
        if cond>=0 && cond<=1
            f(j)=phi(j-1)+(phi(j+1)-phi(j-1))*(cond)*(-4*cond^4+10*cond^3-8*cond^2+cond+2);
        else
            f(j)=phi(j);
        end
    end
    f(1)=phi(1);
    f(N)=phi(N);
end


end

