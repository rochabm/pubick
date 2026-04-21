% crea archivoc para despues dibujarlos
clear 
close all
fileid=fopen('quick_05.txt','w'); 

inicio=0.1;
fin=4;
v=inicio:0.1:fin;

for tf=v
    [numer,exact] = burger_tiempo(256,tf);
    e1=norm(exact-numer,1)/norm(exact,1);
    fprintf(fileid,'%.3f %.3e\n',tf,e1);
end

fclose(fileid);