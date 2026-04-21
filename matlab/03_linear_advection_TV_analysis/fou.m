function [f] = fou(phi)
N=length(phi);
f=zeros(N);

for j=1:N
    f(j)=phi(j);
end

end