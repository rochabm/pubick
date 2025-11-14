function [f] = smart(phi)

N=length(phi);
f=zeros(N);

for j=2:N-1
    if (phi(j+1)-phi(j-1))==0
        f(j)=phi(j);
    else
        cond=(phi(j)-phi(j-1))/(phi(j+1)-phi(j-1));
            
        if cond>=0 && cond<1/6
            f(j)=3*phi(j)-2*phi(j-1);
        else
            if cond>=1/6 && cond<=5/6
                f(j)=(1/8)*(3*phi(j+1)+6*phi(j)-phi(j-1));
            else
                if cond>5/6 && cond<=1
                    f(j)=phi(j+1);
                else
                    f(j)=phi(j);
                end
            end
        end
    end
    f(1)=phi(1);
    f(N)=phi(N);
end   

end