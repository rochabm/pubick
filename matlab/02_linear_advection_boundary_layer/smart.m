function [f_new] = smart(phi,dx,dt,Re)
N=length(phi);
f=zeros(N,1);
f_new=zeros(N,1);

for j=2:N-1

    if (phi(j+1)-phi(j-1))==0
        f(j) = phi(j);
    else
        cond=(phi(j)-phi(j-1))/(phi(j+1)-phi(j-1));
            
        if cond>=0 && cond<1/6
            f(j) = 3*phi(j)-2*phi(j-1);             
        else
            if cond>=1/6 && cond<=5/6
                f(j) = (1/8)*(3*phi(j+1)+6*phi(j)-phi(j-1));                 
            else
                if cond>5/6 && cond<=1
                    f(j) = phi(j+1);                  
                else
                    f(j) = phi(j);
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    conv=(dt/dx)*(f(j)-f(j-1));
    f_new(j) = phi(j) - conv + ((dt)/(Re*dx^2))*(phi(j+1)-2*phi(j)+phi(j-1));

    f_new(1)=phi(1);
    f_new(N)=phi(N);
end

end