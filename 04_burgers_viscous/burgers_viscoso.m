%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Viscous Burgers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all

% ---------------------------------------------

tf  = 10.0;      % tempo final
x0  = 0;         % ponto inicial
xN  = 1;         % ponto final
N   = 256;       % divisao da malha
cfl = 0.001;     % numero de courant(cfl)
dx  = (xN-x0)/N; % espacamento
Re  = 10;        % numero de Reynolds (Re=1/v)
v   = 1/Re;      % viscosidade

N=N+1;

x=zeros(N,1);
%phi_next=zeros(N,1);

dt = cfl*dx;
k = tf/dt;

for i=1:N
    x(i)=x0+(i-1)*dx;
end

% initial condition
phi_in=cond_inicial(x);
phi=phi_in;

% boundary condition
phi(1) = 0;
phi(N) = 0;

% ---------------------------------------------
% exact solution
nterms=100;
phi_exacta=burgers_exacta(x,nterms,v,tf);
% ---------------------------------------------

% ---------------------------------------------
% numerical solution
% ---------------------------------------------

% methods = {@topus, @fou, @pubick, @pubick_calc1, @pubick_calc2, ...
%            @cubick, @sobus, @hpus, @fdhpus, @fdhpus_calc1, ...
%            @smart, @quick};
% 
% method_names = {'topus', 'fou', 'pubick', 'pubick_calc1', 'pubick_calc2', ...
%                 'cubick', 'sobus', 'hpus', 'fdhpus', 'fdhpus_calc1', ...
%                 'smart', 'quick'};

methods = {@fou, @smart, @quick, @topus, @cubick, ...
           @pubick_calc1, @pubick_calc2};

method_names = {'fou', 'smart', 'quick', 'topus', 'cubick', ...
                 'pubick_calc1','pubick_calc2'};
             
% loop over all methods
for m = 1:length(methods)
    fprintf('\nrunning method: %s\n', method_names{m});

    % Reset phi for each run
    phi = phi_in;

    tic;
    % numerical solution
    for j = 1:round(k)
        f = methods{m}(phi, dx, dt, Re);
        % update
        phi = f;
    end
    elapsed = toc;
    fprintf('elapsed time for %s: %.5f seconds\n', method_names{m}, elapsed);
    
    % errors
    numer = phi;
    exact = phi_exacta;
    e1=norm(exact-numer,1)/norm(exact,1);
    e2=norm(exact-numer)/norm(exact);
    einf=norm(exact-numer,Inf)/norm(exact,Inf);
    fprintf('errors: e1 %.8f e2 %.8f ei %.8f\n', e1, e2, einf);
end


%hold on
%plot(x,phi,'b')
%plot(x,phi_exacta,'r')

