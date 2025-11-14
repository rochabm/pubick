%%% CALCULO Y GRAFICO DE BURGER INVISCIDO BESSEL

clear 
close all

% -------------------------------------------------------------------------
% parameters
% -------------------------------------------------------------------------

tf  = 1;         % tiempo final
x0  = 0;         % punto inicial
xN  = 2*pi;      % punto final
N   = 320; %640;       % malha
cfl = 0.001;       % numero de courant(cfl)
dx  = (xN-x0)/N; % espacamento     
v   = 0;         % viscosidad

N=N+1;
x=zeros(N,1);

%phi_next=zeros(N,1);

dt = cfl*dx;
disp(dt);

dt = 1e-5;
cfl = dt/dx;

k = tf/dt;

for i=1:N
    x(i)=x0+(i-1)*dx;
end

% initial condition
phi_in = cond_inicial(x);
phi = phi_in;

% boundary condition
phi(1) = 0;
phi(N) = 0;

% -------------------------------------------------------------------------
% exact solution
% -------------------------------------------------------------------------
fprintf("computing exact solution\n");
n = 200; %500; 
phi_exacta = exacta_bessel(x,n,tf);
% -------------------------------------------------------------------------

% 
% %f1=phi;
% %f2=phi;
% %f3=phi;
% for j=1:round(k) % loop de tiempo
%     %-------------------------------------
%     % Desbloquear el esquema a usar:
%     %-------------------------------------
%     %fa=cubick(f1,dx,dt);
%     %fb=cubick2(f2,dx,dt);
%     f=pubick(phi,dx,dt);
%     %fa=pubick(f1,dx,dt);
%     %fb=pubick2(f2,dx,dt);
%     %f=sobus(phi,dx,dt);
%     %fa=fdhpus(f1,dx,dt);
%     %f=hpus(phi,dx,dt);
%     %f=adbquickest(phi,dx,dt,cfl);
%     %f=smart(phi,dx,dt);
%     %f=topus(phi,dx,dt);
%     %f=fou(phi,dx,dt);
%     %f=epus(phi,dx,dt);
%     %-------------------------------------
%     
%     %%%% actualiza valores;
%     phi=f;
%     %f1=fa;
%     %f2=fb;
%     %f3=fc;
%           
% end


%-------------------------------------------------------
% Define methods and names
%-------------------------------------------------------
methods = {@fou, @smart, @topus, @cubick, @pubick, @quick};
method_names = {'fou', 'smart', 'topus', 'cubick', 'pubick', 'quick',};

fprintf("\ncomputing numerical solutions");

%-------------------------------------------------------
% Loop over methods
%-------------------------------------------------------
for i = 1:length(methods)
    
    tic;
    
    func = methods{i};
    name = method_names{i};
    fprintf("\nrunning method: %s\n",name);
    
    phi = phi_in;

    % Time loop
    for j = 1:round(k)
        % Some schemes (like adbquickest) may need cfl as extra arg
        % Adjust automatically depending on nargin
        if nargin(func) == 4
            f = func(phi, dx, dt, cfl);
        else
            f = func(phi, dx, dt);
        end
        phi = f;
    end

    %-------------------------------------------------------
    % Save final plot silently
    %-------------------------------------------------------
    fig = figure('Visible', 'off');
    
    plot(x, phi_exacta,'r','LineWidth', 4);
    hold on;
    plot(x, phi, ':*b', 'LineWidth',1.5, 'MarkerSize', 4);
    
    title(['Method: ', upper(name)], 'Interpreter', 'none');
    xlabel('x'); ylabel('\phi');
    grid on;
    fname = sprintf('fig_%s.png', name);
    saveas(fig, fname);
    close(fig);

    fprintf('saved: %s\n', fname);
    
    % measure elapsed time
    elapsed = toc;
    fprintf('elapsed time for %s: %.4f seconds\n', method_names{i}, elapsed);
end

fprintf('✅\nall methods executed and plots saved.\n');

% hold on
% set(gca, 'FontSize', 21)
% 
% plot(x,phi_exacta,'r','LineWidth',2)
% plot(x,f,':*b','LineWidth',1.5,'MarkerSize',4)
% plot(x,f1,':*b','LineWidth',1.5,'MarkerSize',4)
% plot(x,f2,':*m','LineWidth',1.5,'MarkerSize',4)
% 
% legend('EXATA','TOPUS','pubick2')
% xlabel('x') 
% ylabel('u')
