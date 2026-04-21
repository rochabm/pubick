function [x,phi_exata,phi_next] = loop_tiempo(tf,N,x0,xN,cfl)

dx=(xN-x0)/N;

N=N+1;
x=zeros(N,1);
phi_next=zeros(N,1);

dt = cfl*dx;
k = tf/dt;

for i=1:N
    x(i)=x0+(i-1)*dx;
end

%---- condicion inicial ----%
%phi_in=problema2(x,N);
phi_in=problema3(x,N);

%%%%%solucion exacta%%%%%%
xt=x-tf; 
%phi_exata=problema2(xt,N);
phi_exata=problema3(xt,N);


phi=phi_in;

%---- solucion numerica ----%
for i=1:round(k) % loop de tiempo

    %Desbloquear esquema a usar:
    %--------------------------------
    %f=topus(phi);
    %f=pubick(phi); %cambiar parametros en pubick.m
    %f=sobus(phi);
    f=hpus(phi); 
    %f=edhpus(phi); %cambiar parametros en edhpus.m
    %f=fdhpus(phi); %cambiar parametros en fdhpus.m
    %f=smart(phi);
    %f=cubick(phi); %cambiar parametros en cubick.m
    %f=adbquickest(phi,cfl);
    %f=fou(phi);
    %------------------------------------
    
    for j=2:N
        phi_next(j)=phi(j)+cfl*(f(j-1)-f(j));
    end
    
    %%% actualiza valores %%%
    phi=phi_next;
        
end


end
