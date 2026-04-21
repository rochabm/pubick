function [phi_in] = problema1(x,N)


phi_in=zeros(N+1,1);

for j=1:N
    phi_in(j) = (sin(pi*x(j)))^4;
end

end