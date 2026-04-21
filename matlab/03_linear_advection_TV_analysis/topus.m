function [f] = topus(phi)
N=length(phi);
f=zeros(N);
alpha=2;
for j=2:N-1
    if (phi(j+1)-phi(j-1))==0
        f(j)=phi(j);
    else
        cond=(phi(j)-phi(j-1))/(phi(j+1)-phi(j-1));
            
        if cond>=0 && cond<=1
            f(j)=phi(j-1)+(phi(j+1)-phi(j-1))*(alpha*cond^4+(-2*alpha+1)*cond^3+((5*alpha-10)/4)*cond^2+((-alpha+10)/4)*cond);
        else
            f(j)=phi(j);
        end
    end
    f(1)=phi(1);
    f(N)=phi(N);
end

end

