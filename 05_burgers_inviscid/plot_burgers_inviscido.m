%%%%%  PLOT BURGERS INVISCIDO %%%%%

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all

%%%%%%%%%%%%%%%%% parametros %%%%%%%%%%%%%%%%%%
tf  = 1;    %tiempo final
x0  = 0;    %punto inicial
xN  = 2*pi;      %punto final
N=600;        %division de la malla
cfl = 0.3;  %numero de courant(cfl)
dx=(xN-x0)/N;      
%v=0;       %viscosidad
%----------------------------------------------
n=500;   %para solucion exacta(Bessel)
%----------------------------------------------

%----------------------------------------------
[An,Ae] = burger_inviscido_func(1,N,cfl,n);
[Bn,Be] = burger_inviscido_func(2,N,cfl,n);
[Cn,Ce] = burger_inviscido_func(3,N,cfl,n);
[Dn,De] = burger_inviscido_func(4,N,cfl,n);
[En,Ee] = burger_inviscido_func(5,N,cfl,n);
[Fn,Fe] = burger_inviscido_func(6,N,cfl,n);
[Gn,Ge] = burger_inviscido_func(7,N,cfl,n);
[Hn,He] = burger_inviscido_func(8,N,cfl,n);
[In,Ie] = burger_inviscido_func(9,N,cfl,n);
[Jn,Je] = burger_inviscido_func(10,N,cfl,n);
[Kn,Ke] = burger_inviscido_func(11,N,cfl,n);
%----------------------------------------------

%----------------------------------------------
N=N+1;
x=zeros(N,1);
%dt = cfl*dx;
%k = tf/dt;
for i=1:N
    x(i)=x0+(i-1)*dx;
end
%%%%%%%%%%%%% solucion exacta %%%%%%%%%%%%%%%%%
phi_exacta=exacta_bessel(x,n,tf);
%----------------------------------------------

hold on
set(gca, 'FontSize', 21)

plot(x,phi_exacta,'-r','LineWidth',2.2)

plot(x,An,':x','LineWidth',2,'MarkerSize',8)
plot(x,Bn,':x','LineWidth',2,'MarkerSize',8)
plot(x,Cn,':x','LineWidth',2,'MarkerSize',8)
plot(x,Dn,':x','LineWidth',2,'MarkerSize',8)
plot(x,En,':x','LineWidth',2,'MarkerSize',8)
plot(x,Fn,':x','LineWidth',2,'MarkerSize',8)
plot(x,Gn,':x','LineWidth',2,'MarkerSize',8)
plot(x,Hn,':x','LineWidth',2,'MarkerSize',8)
plot(x,In,':x','LineWidth',2,'MarkerSize',8)
plot(x,Jn,':x','LineWidth',2,'MarkerSize',8)
plot(x,Kn,':x','LineWidth',2,'MarkerSize',8)

legend('Exata','CUBICK (0,5;0,75)','CUBICK (0,25;0,45)','PUBICK (\mu_1=3/10;\mu_2=5/6)', ...
    'PUBICK (\mu_1=0,493;\mu_2=0,57)','SOBUS','FDHPUS(\theta_1=1,5;\theta_2=0)','HPUS', ...
    'ADBQUICKEST','TOPUS','SMART','FOU')

xlabel('x') 
ylabel('u')
