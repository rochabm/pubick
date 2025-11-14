%Grafica el problema de advección lineal. Usado en 7.1.3 tesis. 

clear 
close all
tf  = 4;
x0  = -1.0;
xN  = 5.0;
cfl=0.5;
 
hold on
% Cambiar esquema en "loop_tiempo.m"
[x50,phi_exata50,phi_next50] = loop_tiempo(tf,50,x0,xN,cfl);
[x100,phi_exata100,phi_next100] = loop_tiempo(tf,100,x0,xN,cfl);
[x200,phi_exata200,phi_next200] = loop_tiempo(tf,200,x0,xN,cfl);

plot(x200,phi_exata200,'-r','LineWidth',1.5) %tv exacta
plot(x50,phi_next50,'--*m','LineWidth',1,'MarkerSize',8) %N=50
plot(x100,phi_next100,'--ob','LineWidth',1,'MarkerSize',8) %N=100
plot(x200,phi_next200,'--diamond g','LineWidth',1,'MarkerSize',8) %N=200

legend('EXATA','N=50','N=100','N=200')
