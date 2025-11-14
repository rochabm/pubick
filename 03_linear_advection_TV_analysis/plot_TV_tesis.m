%Grafica el problema de advección lineal. Usado en 7.1.3 tesis.
%%% Cambiar esquema en "loop_tiempo" %%%

clear 
close all
tiem=zeros(100,1);
TV=zeros(100,1);
TV_exata=zeros(100,1); %tv exata (borrar si no usar)
x0  = -1.0; % dominio: [x0,xN]
xN  = 5.0;
cfl = 0.5;

l=0;
for N=[50 100 200] %aumentar o cambiar N 
    l=l+1;
    for i=1:100
        tempo=i*0.04;
        % Cambiar esquema en "loop_tiempo.m"
        [x,phi_exata,phi] = loop_tiempo(tempo,N,x0,xN,cfl);
    
        tv=0;
        tv_exata=0; % tv exacta
        for j=1:N
            tv=abs(phi(j+1)-phi(j))+tv;
            tv_exata=abs(phi_exata(j+1)-phi_exata(j))+tv_exata; %tv exacta
        end
    
        tiem(i)=tempo;
        TV(i,l)=tv;
        TV_exata(i,l)=tv_exata; %total variation exacta
        
    end
end

hold on

plot(tiem,TV_exata(:,l),'-r','LineWidth',1.5) %tv exacta
plot(tiem,TV(:,1),':*m','LineWidth',1,'MarkerSize',8) %N=50
plot(tiem,TV(:,2),':ob','LineWidth',1,'MarkerSize',8) %N=100
plot(tiem,TV(:,3),':diamond g','LineWidth',1,'MarkerSize',8) %N=200

legend('Exata','N=50','N=100','N=200');
%for i=1:l
%    plot(tiem,TV(:,i),':*')
%end
