function [x,u]=exata2(tf)

N=64;
L=1;
dx=L/N;
N=N+1;
u=zeros(N,1);
x=zeros(N,1);
for j=1:N
    x(j)=(j-1)*dx;
    if x(j)==0
        u(j)=1; 
    else
        xt = x(j)/tf;    
        r2 = sqrt(2.0);
        if (xt < 0.5*(1+r2))
            term = ((-2*xt+sqrt(4*xt+1)-1)/(xt)) + 1;
            u(j) = 0.5*(sqrt(term)+1);      
        else
            u(j)=0.0;
        end
    end
end

end