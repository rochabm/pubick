function [f] = fou(phi)
N=length(phi);
f=zeros(N,1);

for j=1:N
    f(j)=phi(j);
end
f(1)=f(N-1);
f(N)=f(2);
end