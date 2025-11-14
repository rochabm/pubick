function [u] = exacta_bessel(x,n,t)
%v:viscosidad 

m=length(x);
u=zeros(m,1);

for j=1:m
    s=0; 
    for i=1:n
        z=i*t;
        jn=integral_bessel(-z,i);
        s=s+ (jn/z)*sin(i*x(j));  
    end

    u(j)=-2*s;
end

end