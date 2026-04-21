clear
close all

%-------problema 1-------------------
x0=-1; xN=1; N=320; cfl=0.001; tf=1.0; % tf=0.5
%-------problema 2-------------------
%x0=-1; xN=1; N=400; cfl=0.5; tf=0.125;
%-------problema 3-------------------
% x0=-1; xN=5; N=400; cfl=0.5; tf=4;
%------------------------------------

dx=(xN-x0)/N;
N=N+1;

x=zeros(N,1);
phi_next=zeros(N+1,1);

dt = cfl*dx;

dt = 1e-5;
cfl = dt/dx;

k = tf/dt;

for i=1:N
    x(i)=x0+(i-1)*dx;
end

%----------------------------------------------
% condicao inicial
%----------------------------------------------
phi_in=problema1(x,N);
%phi_in=problema2(x,N);
%phi_in=problema3(x,N);

%----------------------------------------------
% exact solution
%----------------------------------------------
xt=x-tf;
phi_exata=problema1(xt,N);
%phi_exata=problema2(xt,N);
%phi_exata=problema3(xt,N);

%----------------------------------------------
% main loop to solve the problem with all methods
%----------------------------------------------

% list of methods to test
% methods = {@topus, @pubick, @sobus, @hpus, @fdhpus, @edhpus, ...
%            @smart, @cubick, @fou};
% 
% method_names = {'topus', 'pubick', 'sobus', 'hpus', 'fdhpus', 'edhpus', ...
%                 'smart', 'cubick', 'fou'};
            
methods = {@fou, @smart, @pubick, @topus, @cubick, @quick};

method_names = {'fou', 'smart', 'pubick', 'topus', 'cubick', 'quick'};

for m = 1:length(methods)
    fprintf('\nrunning method: %s\n', method_names{m});
    
    % Reset phi to initial condition for each run
    phi = phi_in;   
    
    tic;

    % numerical solution time loop
    for i = 1:round(k)
        % select method
        f = methods{m}(phi);  

        % compute solution
        for j = 2:N
            phi_next(j) = phi(j) + cfl * (f(j-1) - f(j));
        end
        
        % periodic condition
        phi_next(1) = phi_next(N);
        
        % update
        phi = phi_next;
    end

    % measure elapsed time
    elapsed = toc;
    fprintf('elapsed time for %s: %.4f seconds\n', method_names{m}, elapsed);
    
    % plot and save silently 
    phi_next_copy = phi_next(1:N);
    phi_exata_copy = phi_exata(1:N);

    fig = figure('Visible', 'off'); % <-- do not show window
    hold on
    plot(x, phi_next_copy, ':*b')
    plot(x, phi_exata_copy, 'r')
    legend('\color{blue}numerical', '\color{red}exact', 'Location', 'Best')
    xlabel('x','FontSize',18)
    ylabel('u','FontSize',18)
    title(sprintf('method: %s', method_names{m}), 'FontSize', 12)
    hold off

    % Save as PNG
    filename = sprintf('fig_%s.png', method_names{m});
    saveas(fig, filename);
    close(fig);
end




