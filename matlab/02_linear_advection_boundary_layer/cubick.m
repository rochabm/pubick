function [f_new] = cubick(phi,dx,dt,Re)
N=length(phi);
f=zeros(N,1);
f_new=zeros(N,1);

%------------- parametros para cubick ---------------------%
a=0.25; b=0.45; c=0.256082165462314; d=0.735437227905742;
%a=0.5; b=0.75; c=0.247622138278052; d=0.674286640691196;
%a=1/3; b=2/3; c=2/3; d=1;
%----------------------------------------------------------
A=(3*c-6*a)/(1+3*a-3*c);
B=(3*a)/(1+3*a-3*c);
C=(-1)/(1+3*a-3*c);
p=(3*B-A^2)/3;
q1=(2*A^3-9*A*B)/27;
%----------------------------------------------------------

for j=2:N-1
    
    if (phi(j+1)-phi(j-1))==0
        f(j)=phi(j);
    else
        cond=(phi(j)-phi(j-1))/(phi(j+1)-phi(j-1));

        if cond>=0 && cond<=1
            q=q1+C*cond;
            auxq=sqrt((q^2)/4+(p^3)/27);
            t=(-q/2+auxq)^(1/3)-(q/2+auxq)^(1/3)-A/3;
            auxf=3*b*(t)*(1-t)^2+3*d*(t^2)*(1-t)+t^3;
            f(j)=phi(j-1)+(phi(j+1)-phi(j-1))*(auxf);
        else
            f(j) = phi(j);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    conv=(dt/dx)*(f(j)-f(j-1));
    f_new(j) = phi(j) - conv + ((dt)/(Re*dx^2))*(phi(j+1)-2*phi(j)+phi(j-1));

    f_new(1)=phi(1);
    f_new(N)=phi(N);
end


end