function [phi,phi_exacta] = burger_inviscido_func(esquema,N,cfl,n)
% BURGER VISCOSO PARA CALCULO DE LA TABLA DE ERRORES (depende de N)

%%%%%%%%%%%%%%%%%% ESQUEMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Escoja opcion de esquema:
%  1: cubick  (parametros (0.5; 0.75) (cambiar parametros en cubick.m))
%  2: cubick2 (parametros (0.25; 0.45) (cambiar parametros en cubick2.m))
%  3: pubick (parametros (3/10; 5/6) (cambiar parametros en pubick.m))
%  4: pubick2 (parametros (0.493; 0.57) (cambiar parametros en pubick2.m))
%  5: sobus
%  6: fdhpus (parametros (1.5; 0) (cambiar parametros en fdhpus.m))
%  7: hpus
%  8: adbquickest
%  9: topus
% 10: smart
% 11: fou
% 12: epus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% parametros %%%%%%%%%%%%%%%%%%%%%%%
tf  = 1;    % tiempo final
x0  = 0;    % punto inicial
xN  = 2*pi; % punto final    
%--------------------------------------------------
%n=500;   % para solucion exacta(Bessel)
%--------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx=(xN-x0)/N;

N=N+1;
x=zeros(N,1);
dt = cfl*dx;
k = tf/dt;
for i=1:N
    x(i)=x0+(i-1)*dx;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% condicion inicial %%%%%%%%%
phi_in=cond_inicial(x);
phi=phi_in;

%%%%%%%%% condicion contorno %%%%%%%%
phi(1)=0;
phi(N)=0;

%%%%%%%%% solucion exacta %%%%%%%%%%%
phi_exacta=exacta_bessel(x,n,tf);


for j=1:round(k) % loop de tiempO
    switch esquema
        case 1
            f=cubick(phi,dx,dt);
        case 2
            f=cubick2(phi,dx,dt);
        case 3
            f=pubick(phi,dx,dt);
        case 4
            f=pubick2(phi,dx,dt);
        case 5
            f=sobus(phi,dx,dt);
        case 6
            f=fdhpus(phi,dx,dt);
        case 7
            f=hpus(phi,dx,dt);
        case 8
            f=adbquickest(phi,dx,dt,cfl);
        case 9
            f=topus(phi,dx,dt);
        case 10
            f=smart(phi,dx,dt);
        case 11
            f=fou(phi,dx,dt);
        case 12
            f=epus(phi,dx,dt);
    end    
    phi=f;  % actualiza valores       
end

%%%%%%% Grafico %%%%%%%%
%hold on
%plot(x,phi,':ob')
%plot(x,phi_exacta,'r')
%-----------------------

end