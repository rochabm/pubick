function [phi,phi_exacta] = adv_linear_suave_func(N)
% ADV. LINEAL SUAVE PARA CALCULO DE LA TABLA DE ERRORES (depende de N)

%%%%%%%%%%%%%%%% parametros %%%%%%%%%%%%%%%%%%%%%%
tf  = 0.5;    %tiempo final
x0  = -1;      %punto inicial
xN  = 1;      %punto final    
cfl = 0.001;  %numero de courant(cfl)
dx=(xN-x0)/N;
%-------------------------------------------------

N=N+1;

x=zeros(N,1);
phi_next=zeros(N+1,1);

dt = cfl*dx;
k = tf/dt;

for i=1:N
    x(i)=x0+(i-1)*dx;
end

%%%%% condicion inicial %%%% (escoger problema)
%----------------------------------------------
phi_in=problema1(x,N);
%----------------------------------------------

%%%%%% solucion exacta %%%%% (escoger problema)
%----------------------------------------------
xt=x-tf;
phi_exacta=problema1(xt,N);
%----------------------------------------------

phi=phi_in;

%%%%% solucion numerica %%%%%
for i=1:round(k) % loop de tiempo
    %phi(1)=phi(N-1);
    %phi(N)=phi(2);
    %-----------------------------------
    %Desbloquear f segun esquema a usar: 
    %-----------------------------------
    %f=topus(phi);
    %f=pubick(phi); %cambiar parametros en pubick.m
    %f=sobus(phi);
    %f=quick(phi);
    %f=hpus(phi); 
    %f=fdhpus(phi); %cambiar parametros en fdhpus.m
    %f=edhpus(phi); %cambiar parametros en edhpus.m
    %f=smart(phi);
    f=cubick(phi); %cambiar parametros en cubick.m 
    %f=adbquickest(phi,cfl);
    %f=fou(phi);
    %f=quick(phi);
    %-----------------------------------
    
    for j=2:N
        phi_next(j)=phi(j)+cfl*(f(j-1)-f(j));
    end
    
    %%% condicion periodica
    phi_next(1)=phi(N);
    %phi_next(N+1)=phi_next(2);
    

    %%% actualiza valores %%%
    phi=phi_next;
end
 

