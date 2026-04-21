function [phi,phi_exacta] = burger_viscoso_func(N)
% BURGER VISCOSO PARA CALCULO DE LA TABLA DE ERRORES (depende de N)

%%%%%%%%%%%%%%%% parametros %%%%%%%%%%%%%%%%%%%%%%
tf  = 0.8;    %tiempo final
x0  = 0;      %punto inicial
xN  = 1;      %punto final    
cfl = 0.001;  %numero de courant(cfl)
dx=(xN-x0)/N;
Re=10;        %numero de Reynolds (Re=1/v)
v=1/Re;       %viscosidad
%-------------------------------------------------
n=100;   %para solucion exacta(series de fourier)
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
phi(N)=0;

%%%%%%%%% solucion exacta %%%%%%%%%%
phi_exacta=burgers_exacta(x,n,v,tf);

for j=1:round(k) % loop de tiempo
    %---------------------------------
    % Desbloquear el esquema a usar:
    %---------------------------------
    %f=topus(phi,dx,dt,Re);
    %f=fou(phi,dx,dt,Re);
    %f=pubick(phi,dx,dt,Re); %modifique parametros en "pubick.m"
    %f=cubick(phi,dx,dt,Re); %modifique parametros en "cubick.m"
    %f=sobus(phi,dx,dt,Re);
    f=quick(phi,dx,dt,Re);
    %f=hpus(phi,dx,dt,Re);
    %f=fdhpus(phi,dx,dt,Re); %modifique parametros en "fdhpus.m"
    %f=adbquickest(phi,dx,dt,Re,cfl);
    %f=smart(phi,dx,dt,Re);
    %---------------------------------
    
    phi=f;         
end

end