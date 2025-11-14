function [f] = adbquickest(phi,CFL)
N=length(phi);
f=zeros(N);

numa = ((1/3)-(1/2)*abs(CFL)+(1/6)*CFL^2);
dena = (2*((7/12)-(1/2)*CFL-(1/4)*abs(CFL)+(1/6)*CFL^2));
a    = numa/dena;

numb = ((1/6)*CFL^2-(2/3)-(1/2)*abs(CFL)+CFL);
denb = (2*((-5/12)-(1/4)*abs(CFL)+(1/6)*CFL^2+(1/2)*CFL));
b    = numb/denb;

for j=2:N-1
    if (phi(j+1)-phi(j-1))==0
        f(j)=phi(j);
    else
        cond=(phi(j)-phi(j-1))/(phi(j+1)-phi(j-1));
         A=phi(j)-phi(j-1);
         B=phi(j+1)-phi(j-1);
         C=phi(j-1);
            
        if cond>=a && cond<=b
            f(j) =C+A+0.5*(1-abs(CFL))*(B-A)-(1/6)*(1-CFL^2)*(B-2*A);
        else
            if cond>=0 && cond<a
                f(j)=C+(2-CFL)*A;
            else
                if cond>b && cond<=1
                    f(j)=C+B-B*CFL+A*CFL;
                else
                    f(j)=phi(j);
                end
            end
        end
    end
    f(1)=phi(1);
    f(N)=phi(N);
end

end


 


