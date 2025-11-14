function [phi_in] = problema3(x,N)


phi_in=zeros(N,1);

for j=1:N
    if x(j)>=-1/3 && x(j)<=1/3
        phi_in(j)=1;
    else
        phi_in(j)=0;
    end
end

end

