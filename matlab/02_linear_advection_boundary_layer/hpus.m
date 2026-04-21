function [f_new] = hpus(phi,dx,dt,Re)
N=length(phi);
f=zeros(N,1);
f_new=zeros(N,1);

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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    conv=(dt/dx)*(f(j)-f(j-1));
    f_new(j) = phi(j) - conv + ((dt)/(Re*dx^2))*(phi(j+1)-2*phi(j)+phi(j-1));

    f_new(1)=phi(1);
    f_new(N)=phi(N);
end

end