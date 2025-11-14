clear all
close all
clc
%malla=200, cfl=0.6

load cubick_N200_cfl06_caso1.txt
load minmod_N200_cfl06_caso1.txt
load pubick_N200_cfl06_caso1.txt
load smart_N200_cfl06_caso1.txt
load superbee_N200_cfl06_caso1.txt
load topus_N200_cfl06_caso1.txt

load Exata_caso1.txt

A=cubick_N200_cfl06_caso1;
B=pubick_N200_cfl06_caso1;
C=superbee_N200_cfl06_caso1;
D=minmod_N200_cfl06_caso1;
E=smart_N200_cfl06_caso1;
F=topus_N200_cfl06_caso1;

G=Exata_caso1;

hold on
set(gca, 'FontSize', 21)

plot(G(:,1),G(:,2),'r','LineWidth',2.2)

plot(A(:,1),A(:,2),':x','LineWidth',2,'MarkerSize',8)
plot(B(:,1),B(:,2),':x','LineWidth',2,'MarkerSize',8)
plot(C(:,1),C(:,2),':x','LineWidth',2,'MarkerSize',8)
plot(D(:,1),D(:,2),':x','LineWidth',2,'MarkerSize',8)
plot(E(:,1),E(:,2),':x','LineWidth',2,'MarkerSize',8)
plot(F(:,1),F(:,2),':x','LineWidth',2,'MarkerSize',8)


legend('Exata','CUBICK (0,5;0,75)','PUBICK (\mu_1=3/10;\mu_2=5/6)', ...
    'Superbee','MINMOD','SMART','TOPUS')

xlabel('x') 
ylabel('\rho')

%ylim([1 4.8]);
