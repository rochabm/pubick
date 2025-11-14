clear all
close all
clc
%malla=300, cfl=0.6, tf=1

load cubick_cfl06_N2500_t5.txt
load minmod_cfl06_N2500_t5.txt
load pubick_cfl06_N2500_t5.txt
load smart_cfl06_N2500_t5.txt
load superbee_cfl06_N2500_t5.txt
load topus_cfl06_N2500_t5.txt

load Reference_Upwind.txt

A=cubick_cfl06_N2500_t5;
B=pubick_cfl06_N2500_t5;
C=superbee_cfl06_N2500_t5;
D=minmod_cfl06_N2500_t5;
E=smart_cfl06_N2500_t5;
F=topus_cfl06_N2500_t5;

G=Reference_Upwind;

hold on
set(gca, 'FontSize', 21)

plot(G(:,1),G(:,2),'r','LineWidth',1)

%plot(A(:,1),A(:,2),'squareb','LineWidth',1.2,'MarkerSize',4); legend('Exata','CUBICK (0,5;0,75)');
%plot(B(:,1),B(:,2),'squareb','LineWidth',1.2,'MarkerSize',4); legend('Exata','PUBICK (\mu_1=3/10;\mu_2=5/6)');
%plot(C(:,1),C(:,2),'squareb','LineWidth',1.2,'MarkerSize',4); legend('Exata','Superbee');
%plot(D(:,1),D(:,2),'squareb','LineWidth',1.2,'MarkerSize',4); legend('Exata','MINMOD');
%plot(E(:,1),E(:,2),'squareb','LineWidth',1.2,'MarkerSize',4); legend('Exata','SMART');
plot(F(:,1),F(:,2),'squareb','LineWidth',1.2,'MarkerSize',4); legend('Exata','TOPUS');

xlabel('x') 
ylabel('\rho')

ylim([0.88 1.7]);
