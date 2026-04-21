function [u] = burgers_exacta(x,n,v,t)
%v:viscosidad 

m=length(x);
u=zeros(m,1);

a0=integral_a0(v);

for j=1:m
    s1=0; %s1: sumatoria numerador
    s2=0; %s2: sumatoria denominador
    for i=1:n
        an=integral_an(v,i);
        s1=s1+an*(exp(-i^2*pi^2*v*t))*i*sin(i*pi*x(j));
        s2=s2+an*(exp(-i^2*pi^2*v*t))*cos(i*pi*x(j));   
    end

    u(j)=(2*pi*v)*((s1)/(a0+s2));
end

end