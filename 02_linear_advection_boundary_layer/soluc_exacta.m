function [u] = camada_limite_exacta(x,v)
%v:viscosidad 

m=length(x);
u=zeros(m,1);

for j=1:m

    u(j)=(1-exp(x(j)/v))/(1-exp(1/v));

end

end