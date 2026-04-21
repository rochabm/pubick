function [f_new] = adbquickest(phi,dx,dt,Re,cfl)
N=length(phi);
f=zeros(N,1);
f_new=zeros(N,1);

%---------------------- parametros -----------------------
alfad  = (1/2)*(1-abs(cfl))-(1/6)*(1-cfl^2);
alfau  = 1-(1/2)*(1-abs(cfl))+(1/3)*(1-cfl^2);
alfar  = (1/6)*(1-cfl^2);

numa = (2-(3*abs(cfl))+(cfl^2));
dena = (7-(6*cfl)-(3*abs(cfl))+(2*cfl^2));
a    = numa/dena;

numb = ((cfl^2)-4-(3*abs(cfl))+(6*cfl));
denb = (-5-(3*abs(cfl))+(2*cfl^2)+(6*cfl));
b    = numb/denb;
%---------------------------------------------------------

for j=2:N-1

    if (phi(j+1)-phi(j-1))==0
        f(j)=phi(j);
    else
        cond=(phi(j)-phi(j-1))/(phi(j+1)-phi(j-1));
               
        if 0<=cond && cond<a
            f(j) = (2-cfl)*phi(j)-(1-cfl)*phi(j-1);
        else
            if a<=cond && cond<=b
                f(j) = alfad*phi(j+1)+alfau*phi(j)-alfar*phi(j-1);
            else
                if b<cond && cond<=1
                    f(j) = (1.0-cfl)*phi(j+1)+cfl*phi(j);
                else
                    f(j) = phi(j);
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    conv=(dt/dx)*(f(j)-f(j-1));
    f_new(j) = phi(j) - conv + ((dt)/(Re*dx^2))*(phi(j+1)-2*phi(j)+phi(j-1));

    f_new(1)=phi(1);
    f_new(N)=phi(N);
end

end