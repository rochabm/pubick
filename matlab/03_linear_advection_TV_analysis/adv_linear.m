clear
close all
%-------problema 2-------------------
x0=-1; xN=1; N=400; cfl=0.5; tf=0.125;
%-------problema 3-------------------
% x0=-1; xN=5; N=400; cfl=0.5; tf=4;
%------------------------------------

dx=(xN-x0)/N;
N=N+1;

x=zeros(N,1);
phi_next=zeros(N,1);

dt = cfl*dx;
k = tf/dt;

for i=1:N
    x(i)=x0+(i-1)*dx;
end

%%%%% condicion inicial %%%% (escoger problema)
%----------------------------------------------
phi_in=problema2(x,N);
%phi_in=problema3(x,N);
%----------------------------------------------

%%%%%% solucion exacta %%%%% (escoger problema)
%----------------------------------------------
xt=x-tf;
phi_exata=problema2(xt,N);
%phi_exata=problema3(xt,N);
%----------------------------------------------

phi=phi_in;

%%%%% solucion numerica %%%%%
for i=1:round(k) % loop de tiempo
    %-----------------------------------
    %Desbloquear f segun esquema a usar: 
    %-----------------------------------
    %f=topus(phi);
    %f=pubick(phi); %cambiar parametros en pubick.m
    %f=sobus(phi);
    %f=hpus(phi); 
    %f=fdhpus(phi); %cambiar parametros en fdhpus.m
    %f=edhpus(phi); %cambiar parametros en edhpus.m
    %f=smart(phi);
    %f=cubick(phi); %cambiar parametros en cubick.m 
    %f=adbquickest(phi,cfl);
    %f=fou(phi);
    f=quick(phi);
    %-----------------------------------
    
    for j=2:N
        phi_next(j)=phi(j)+cfl*(f(j-1)-f(j));
    end
    
    %%% actualiza valores %%%
    phi=phi_next;
        
end

hold on
set(gca, 'FontSize', 21)


plot(x,phi_exata,'r','LineWidth',2)
plot(x,phi_next,':*b','LineWidth',1.5,'MarkerSize',8)
legend('Exact','TOPUS')
xlabel('x') 
ylabel('u')



