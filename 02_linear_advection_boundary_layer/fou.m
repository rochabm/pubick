function [f_new] = fou(phi,dx,dt,Re)
N=length(phi);
f=zeros(N,1);
f_new=zeros(N,1);

for j=2:N-1

    f(j) = phi(j);       

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    conv=(dt/dx)*(f(j)-f(j-1));
    f_new(j) = phi(j) - conv + ((dt)/(Re*dx^2))*(phi(j+1)-2*phi(j)+phi(j-1));

    f_new(1)=phi(1);
    f_new(N)=phi(N);
end

end