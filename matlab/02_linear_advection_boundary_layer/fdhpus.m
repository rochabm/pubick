function [f_new] = fdhpus(phi,dx,dt,Re)
N=length(phi);
f=zeros(N,1);
f_new=zeros(N,1);
%--------- parametros fdhpus ------------
a=1.5; b=0; % a=\theta_1, b=\theta_2 
%a=2; b=0; %(wei)
%----------------------------------------

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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    conv=(dt/dx)*(f(j)-f(j-1));
    f_new(j) = phi(j) - conv + ((dt)/(Re*dx^2))*(phi(j+1)-2*phi(j)+phi(j-1));

    f_new(1)=phi(1);
    f_new(N)=phi(N);
end

end