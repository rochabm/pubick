clear 
close all

tf= 0.5;
x0= 0;
xN= 1;
N=50;
c=1;
cfl = 0.2;
dx=(xN-x0)/N;
dt = cfl*dx; %
%--------------------------
N=N+1;
%u=zeros(N,1);

x=zeros(N,1);
for i=1:N
    x(i)=x0+(i-1)*dx;
end

u=condicion_inicial(x);
Exacta=exata(x,tf);

%A=cubick(u,dt,dx,tf,c);
B=cubick2(u,dt,dx,tf,c);
C=pubick(u,dt,dx,tf,c);
D=pubick2(u,dt,dx,tf,c);
%E=sobus(u,dt,dx,tf,c);
%F=fdhpus(u,dt,dx,tf,c);
%G=hpus(u,dt,dx,tf,c);
%H=adbquickest(u,dt,dx,tf,c,cfl);
I=topus(u,dt,dx,tf,c);
J=smart(u,dt,dx,tf,c);
%K=fou(u,dt,dx,tf,c);
L=quick(u,dt,dx,tf,c);

hold on
set(gca, 'FontSize', 21)

plot(x,Exacta,'r','LineWidth',2.2)

%plot(x,A,':x','LineWidth',2,'MarkerSize',8) %cubick1
plot(x,B,':x','LineWidth',2,'MarkerSize',8) %cubick2
plot(x,C,':x','LineWidth',2,'MarkerSize',8) %pubick1
plot(x,D,':x','LineWidth',2,'MarkerSize',8) %pubick2
%plot(x,E,':x','LineWidth',2,'MarkerSize',8) %sobus
%plot(x,F,':x','LineWidth',2,'MarkerSize',8) %fdhpus
%plot(x,G,':x','LineWidth',2,'MarkerSize',8) %hpus
%plot(x,H,':x','LineWidth',2,'MarkerSize',8) %adbquickest
plot(x,I,':x','LineWidth',2,'MarkerSize',8) %topus
plot(x,J,':x','LineWidth',2,'MarkerSize',8) %smart
%plot(x,K,':x','LineWidth',2,'MarkerSize',8) % fou
plot(x,L,':x','LineWidth',2,'MarkerSize',8) % quick


%legend('Exata','CUBICK (0,5;0,75)','CUBICK (0,25;0,45)','PUBICK (\mu_1=3/10;\mu_2=5/6)', ...
%    'PUBICK (\mu_1=0,493;\mu_2=0,57)','SOBUS','FDHPUS(\theta_1=1,5;\theta_2=0)','HPUS', ...
%    'ADBQUICKEST','TOPUS','SMART','FOU')

legend('Exact','CUBICK (0.25,0.45)','PUBICK (\mu_1=3/10,\mu_2=5/6)', ...
    'PUBICK (\mu_1=0.493,\mu_2=0.57)','TOPUS','SMART','QUICK')

xlabel('x') 
ylabel('u')

ylim([-0.05 1.05]);







