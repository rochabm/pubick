%%% CALCULO, GRAFICO Y TIEMPO DE SIMULACION DE BURGER VISCOSO 

clear 
close all

%%%%%%%%%%%%%%%%%%%%
tf  = 0.5;    %tiempo final
x0  = 0;      %punto inicial
xN  = 1;      %punto final
N=512;        %division de la malla
cfl = 0.01;  %nmero de courant(cfl)
dx=(xN-x0)/N;
Re=100;        %nmero de Reynolds (Re=1/v)
v=1/Re;       %viscosidad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=N+1;

x=zeros(N,1);
%phi_next=zeros(N,1);

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

tic;
for j=1:round(k) % loop de tiempo
    %-------------------------------------
    % Desbloquear el esquema a usar:
    %-------------------------------------
    %f=topus(phi,dx,dt,Re);
    %f=fou(phi,dx,dt,Re);
    %f=pubick(phi,dx,dt,Re);
    %f=cubick(phi,dx,dt,Re);
    %f=sobus(phi,dx,dt,Re);
    %f=hpus(phi,dx,dt,Re);
    %f=fdhpus(phi,dx,dt,Re);
    %f=adbquickest(phi,dx,dt,Re,cfl);
    f=smart(phi,dx,dt,Re);
    %-------------------------------------
    
    %%%% actualiza valores;
    phi=f;
          
end

elapsedTime = toc;
fprintf('El tiempo de simulación es: %.5f segundos.\n', elapsedTime);


hold on

plot(x,phi,'b')
plot(x,phi_exacta,'r')

