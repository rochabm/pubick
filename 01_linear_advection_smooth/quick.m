function [f] = quick(phi)
N=length(phi);
f=zeros(N,1);

for j=2:N-1
    if (phi(j+1)-phi(j-1))==0
        f(j)=phi(j);
    else
        cond=(phi(j)-phi(j-1))/(phi(j+1)-phi(j-1));
            
        if cond>=0 && cond<=1
            f(j)=phi(j-1)+(phi(j+1)-phi(j-1))*((6/8)*cond+3/8);
        else
            f(j)=phi(j);
        end
    end
end 
f(1)=f(N-1);
f(N)=f(2);
end