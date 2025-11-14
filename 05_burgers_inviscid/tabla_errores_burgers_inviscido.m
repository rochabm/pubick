%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TABLA DE ERRORES DE CONVERGENCIA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% numer: solucion numerica
% exact: solucion exacta
% e1   : error 1 (norma 1)
% e2   : error 2 (norma 2)
% einf : error inf (norma infinita)

esquema=12;

l=0;
for N=[50 100 200 400] %aumentar N 
    l=l+1;

    [numer,exact]=burger_inviscido_func(esquema,N);

    e1=norm(exact-numer,1)/norm(exact,1);
    e2=norm(exact-numer)/norm(exact);
    einf=norm(exact-numer,Inf)/norm(exact,Inf);
    q=[e1 e2 einf];
    if l==1
        p=zeros(1,3);
    else
        p=log(q./q_ant)/log(1/2); % calculo del orden
    end

    Salida=[e1 p(1) e2 p(2) einf p(3)];

    fprintf('%d\t|  %.3e  |  %.3f  |  %.3e  |  %.3f  |  %.3e  |  %.3f  |\n',N,Salida);
    q_ant=q;


end
