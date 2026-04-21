clear all
close all
clc
%malla=50, cfl=0.6, tf=0.2

load cubick_N50_cfl06_t02.txt
load minmod_N50_cfl06_t02.txt
load pubick_N50_cfl06_t02.txt
load smart_N50_cfl06_t02.txt
load superbee_N50_cfl06_t02.txt
load topus_N50_cfl06_t02.txt

load exata_cfl06.txt

A=cubick_N50_cfl06_t02;
B=pubick_N50_cfl06_t02;
C=superbee_N50_cfl06_t02;
D=minmod_N50_cfl06_t02;
E=smart_N50_cfl06_t02;
F=topus_N50_cfl06_t02;

G=exata_cfl06;

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

ylim([0.1 1.08]);
