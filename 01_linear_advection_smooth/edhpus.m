function [f] = edhpus(phi)
%Esquema Hermite grado 8

N=length(phi);
f=zeros(N,1);
a=1; %theta 1 (a=2 wei)
b=1; %theta 2 (b=0 wei)
c=0;
d=(95/8)-9;
e=0;
for j=2:N-1
    if (phi(j+1)-phi(j-1))==0
        f(j)=phi(j);
    else
        cond=(phi(j)-phi(j-1))/(phi(j+1)-phi(j-1));
            
        if cond>=0 && cond<=1
            %f(j)=phi(j-1)+(phi(j+1)-phi(j-1))*(cond)*(-4*cond^4+10*cond^3-8*cond^2+cond+2);
            aux=(72*a-72*b+4*c-32*d+4*e-192)*cond^8+...
                (-316*a+260*b-18*c+128*d-14*e+840)*cond^7+...
                (558*a-362*b+33*c-200*d+19*e-1468)*cond^6+...
                (-501*a+243*b-31.5*c+152*d-12.5*e+1290)*cond^5+...
                (234*a-79*b+16.5*c-56*d+4*e-579)*cond^4+...
                (-48*a+10*b-4.5*c+8*d-0.5*e+110)*cond^3+...
                (0.5*c)*cond^2+...
                (a)*cond;
            f(j)=phi(j-1)+(phi(j+1)-phi(j-1))*(aux);
        else
            f(j)=phi(j);
        end
    end
end
f(1)=f(N-1);
f(N)=f(2);

end