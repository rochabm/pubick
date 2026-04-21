function [u] = condicion_inicial(x)

N=length(x);
u=zeros(N,1);
for j=1:N
    if (x(j)==0)
       u(j)=1;
    else
        u(j)=0;
    end
end

end