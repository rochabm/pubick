%%% Cambiar esquema en "loop_tiempo" %%%

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

clear 
close all
tiem=zeros(100,1);
TV=zeros(100,1);
TV_exata=zeros(100,1); %tv exata (borrar si no usar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0  = 0;    % punto inicial
xN  = 2*pi; % punto final    
cfl = 0.3;  % numero de courant(cfl)
%--------------------------------------------------
%n=500;   % para solucion exacta(Bessel)
%--------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l=0;
for N=[100 200 400] %aumentar N 
    l=l+1;

    dx=(xN-x0)/N;
    x=zeros(N+1,1);
    dt = cfl*dx;
    
    
    for j=1:(N+1)
        x(j)=x0+(j-1)*dx;
    end
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    phi_exata=exacta_bessel(x,200,0.25);  
         
    %%%%%%%%% condicion inicial %%%%%%%%%    
    phi_in=cond_inicial(x);
    phi=phi_in;

    %%%%%%%%% condicion contorno %%%%%%%%
    phi(1)=0;
    phi(N+1)=0;    

    for i=1:100

        %%%%%%%%%%%%%%%% parametros %%%%%%%%%%%%%%%%%%%%%%%
        tf  = i*0.0025;    % tiempo final 
        k = tf/dt;
        %--------------------------------------------------

        for j=1:round(k) % loop de tiempO
            f=topus(phi,dx,dt);
            phi=f;  % actualiza valores       
        end
    
        tv=0;
        tv_exata=0; % tv exacta
        for j=1:N
            tv=abs(phi(j+1)-phi(j))+tv;
            tv_exata=abs(phi_exata(j+1)-phi_exata(j))+tv_exata; %tv exacta
        end
    
        tiem(i)=tf;
        TV(i,l)=tv;
        TV_exata(i,l)=tv_exata; %total variation exacta
        
    end
end




hold on

plot(tiem,TV_exata(:,l),'-r','LineWidth',1.5) %tv exacta
plot(tiem,TV(:,3),':*m','LineWidth',1,'MarkerSize',8) %N=400
plot(tiem,TV(:,1),':ob','LineWidth',1,'MarkerSize',8) %N=100
plot(tiem,TV(:,2),':diamond g','LineWidth',1,'MarkerSize',8) %N=200

legend('Exata','N=400','N=100','N=200');
%for i=1:l
%    plot(tiem,TV(:,i),':*')
%end