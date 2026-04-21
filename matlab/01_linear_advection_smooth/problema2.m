function [phi_in] = problema2(x,N)
%7.1.2 Condição Inicial na Forma de um “W”

phi_in=zeros(N,1);

for j=1:N
    if x(j)>=-1.0 && x(j)<= 0.0
        phi_in(j)=0;
    else
        if x(j) >= 0.0 && x(j) <= 0.2
            phi_in(j)=1;
        else
            if x(j) > 0.2 && x(j) <= 0.4
                phi_in(j)=4*x(j)-3/5;
            else
                if x(j) > 0.4 && x(j) <= 0.6
                    phi_in(j)=-4*x(j)+13/5;
                else
                    if x(j) > 0.6 && x(j) <= 0.8
                        phi_in(j)=1;
                    else
                        if x(j) > 0.8 && x(j) <= 1
                            phi_in(j)=0;
                        end
                    end
                end
            end
        end
    end
end


end

