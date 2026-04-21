function [phi,phi_exacta] = camada_limite_func(N)
% BURGER VISCOSO PARA CALCULO DE LA TABLA DE ERRORES (depende de N)

%%%%%%%%%%%%%%%% parametros %%%%%%%%%%%%%%%%%%%%%%
tf  = 0.5;    %tiempo final
x0  = 0;      %punto inicial
xN  = 1;      %punto final    
cfl = 0.01;  %numero de courant(cfl)
dx=(xN-x0)/N;
Re=100;        %numero de Reynolds (Re=1/v)
v=1/Re;       %viscosidad
%-------------------------------------------------

N=N+1;

x=zeros(N,1);

dt = cfl*dx;
k = tf/dt;

for i=1:N
    x(i)=x0+(i-1)*dx;
end

%%%%%%%% condicion inicial %%%%%%%%
phi_in=cond_inicial(x);
phi=phi_in;

%%%%%%%% condicion contorno %%%%%%%%
phi(1)=0;
phi(N)=1;

%%%%%%%%% solucion exacta %%%%%%%%%%
phi_exacta=soluc_exacta(x,v);

for j=1:round(k) % loop de tiempo
    %---------------------------------
    % Desbloquear el esquema a usar:
    %---------------------------------
    %f=topus(phi,dx,dt,Re);
    %f=fou(phi,dx,dt,Re);
    %f=pubick(phi,dx,dt,Re); %modifique parametros en "pubick.m"
    %f=cubick(phi,dx,dt,Re); %modifique parametros en "cubick.m"
    %f=sobus(phi,dx,dt,Re);
    %f=hpus(phi,dx,dt,Re);
    %f=fdhpus(phi,dx,dt,Re); %modifique parametros en "fdhpus.m"
    %f=adbquickest(phi,dx,dt,Re,cfl);
    %f=smart(phi,dx,dt,Re);
    f=quick(phi,dx,dt,Re);
    %---------------------------------
    
    phi=f;         
end
%hold on
%plot(x,phi_exacta);
%plot(x,phi);

end