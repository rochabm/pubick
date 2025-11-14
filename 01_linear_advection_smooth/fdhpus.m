function [f] = fdhpus(phi)

N=length(phi);
f=zeros(N,1);
a=1.5; %theta 1 (a=2 wei)
b=0; %theta 2 (b=0 wei)
for j=2:N-1
    if (phi(j+1)-phi(j-1))==0
        f(j)=phi(j);
    else
        cond=(phi(j)-phi(j-1))/(phi(j+1)-phi(j-1));
            
        if cond>=0 && cond<=1
            f(j)=phi(j-1)+(phi(j+1)-phi(j-1))*(4*(a+b-3)*cond^5-2*(6*a+4*b-17)*cond^4+(13*a+5*b-34)*cond^3-(6*a+b-13)*cond^2+a*cond);
        else
            f(j)=phi(j);
        end
    end
end
f(1)=f(N-1);
f(N)=f(2);

end

