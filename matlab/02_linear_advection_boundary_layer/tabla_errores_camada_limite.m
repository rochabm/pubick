%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TABLA DE ERRORES DE CONVERGENCIA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% numer: solucion numerica
% exact: solucion exacta
% e1   : error 1 (norma 1)
% e2   : error 2 (norma 2)
% einf : error inf (norma infinita)


clear 
close all

l=0;
for N=[512 1024 2048 4096] %aumentar N 
    l=l+1;

    % Cambie esquema en "camada_limite_func.m"
    [numer,exact]=camada_limite_func(N);

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
