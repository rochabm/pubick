function [phi_in] = cond_inicial(x)

N=length(x);
phi_in=zeros(N,1);

for j=1:N
    phi_in(j)=sin(pi*(x(j)));   
end

end

