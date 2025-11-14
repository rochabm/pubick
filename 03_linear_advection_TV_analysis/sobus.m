function [f] = sobus(phi)
N=length(phi);
f=zeros(N);

for j=2:N-1
    if (phi(j+1)-phi(j-1))==0
        f(j)=phi(j);
    else
        cond=(phi(j)-phi(j-1))/(phi(j+1)-phi(j-1));
            
        if cond>=0 && cond<=1
            f(j)=phi(j-1)+(phi(j+1)-phi(j-1))*(-(sqrt(3)/2)*cond+(1/2+sqrt(3)/3)*(-(1/2)*(3-sqrt(3))+3*sqrt((1/6)*(2-sqrt(3))+(sqrt(3)/3)*cond)));
        else
            f(j)=phi(j);
        end
    end
    f(1)=phi(1);
    f(N)=phi(N);
end

end
