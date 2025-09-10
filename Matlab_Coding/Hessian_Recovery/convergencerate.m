clc
clear all
for j = 1:7
    if j==1
        disp(['Example' num2str(j) ':the example is a smooth function'])
    else
        disp(['Example ' num2str(j) ': the degree of the polynomial function is ' num2str(j)]);
    end
    M = 4;
    error = zeros(M,1);
    h = zeros(M,1);
    runtime = zeros(M,1);
    Gh_error = zeros(M,1);
    Gin_error = zeros(M,1);
    Hin_error = zeros(M,1);
    Hh_error = zeros(M,1);
    Dofs = zeros(M,1);
    for i=1:M
        N = 2^(i+3);
        h(i) = 1/(2^(i+3));
        tic
        [node,elem] = cartesian(N,N);
        Dofs(i) = size(node,1);
        [Gux,Guy,Hxx,Hxy,Hyx,Hyy] = hessianrecovery2(node,elem);
        x = node(:,1);
        y = node(:,2);
        [u,ux,uy,uxx,uxy,uyx,uyy] = exact_function(x,y,num2str(j));
        Gh_ux = Gux*u;
        Gh_uy = Guy*u;
        errorx =abs(ux-Gh_ux);
        errory= abs(uy-Gh_uy);
        L = 0.1;
        interior = find(node(:,1)<1-L&node(:,1)>L&node(:,2)<1-L&node(:,2)>L);
        in_errorX = max(errorx(interior(:,1),1));
        in_errorY = max(errory(interior(:,1),1));
        Gh_error(i) = max(max(errorx),max(errory));
        Gin_error(i) = max(in_errorX,in_errorY);
        % Hessian recovery
        Hh_uxx = Hxx*u;
        Hh_uxy = Hxy*u;
        Hh_uyx = Hyx*u;
        Hh_uyy = Hyy*u;
        errorxx = abs(uxx-Hh_uxx);
        errorxy = abs(uxy-Hh_uxy);
        erroryx = abs(uyx-Hh_uyx);
        erroryy = abs(uyy-Hh_uyy);
        in_errorXX = max(errorxx(interior(:,1),1));
        in_errorXY = max(errorxy(interior(:,1),1));
        in_errorYX = max(erroryx(interior(:,1),1));
        in_errorYY = max(erroryy(interior(:,1),1));
        Hin_error(i) = max([in_errorXX,in_errorXY,in_errorYX,in_errorYY]);
        Hh_error(i) = max([max(errorxx),max(errorxy),max(erroryx),max(erroryy)]);
        runtime(i) = toc;
    end

    figure(3)
    loglog(Dofs,Gin_error,'r-*',Dofs,Gh_error,'b-*', ...
        Dofs,Hin_error,'m-+',Dofs,Hh_error,'-+');
    legend('$\|G_hu-\nabla u\|_{\infty,\Omega_{1,h}}$','$\|G_hu-\nabla u\|_{\infty,\Omega_h}$', ...
        '$\|H_hu-\Delta u\|_{\infty,\Omega_{1,h}}$','$\|H_hu-\Delta u\|_{\infty,\Omega_h}$','Interpreter','latex','Location','best');
    xlabel('Number of Degree of Freedoms');
    ylabel('Error');
    fig = gcf;
    fig.PaperPositionMode ='auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(['CRateAllOperator_',num2str(j)],'-dpdf');

    figure(4)
    loglog(Dofs,Gin_error,'r-*',Dofs,Gh_error,'b-*', ...
        Dofs,Hin_error,'m-+');
    legend('$\|G_hu-\nabla u\|_{\infty,\Omega_{1,h}}$','$\|G_hu-\nabla u\|_{\infty,\Omega_h}$', ...
        '$\|H_hu-\Delta u\|_{\infty,\Omega_{1,h}}$','Interpreter','latex','Location','best');
    xlabel('Number of Degree of Freedoms');
    ylabel('Error');
    fig = gcf;
    fig.PaperPositionMode ='auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(['CRateOperator_',num2str(j)],'-dpdf');



    figure(5)
    loglog(Dofs,Gin_error,'r-*',Dofs,Gh_error,'b-*');
    legend('$\|G_hu-\nabla u\|_{\infty,\Omega_{1,h}}$','$\|G_hu-\nabla u\|_{\infty,\Omega_h}$', ...
        'Interpreter','latex','Location','best');
    xlabel('Number of Degree of Freedoms');
    ylabel('${\|G_hu-\nabla u\|}$','Interpreter','latex');
    fig = gcf;
    fig.PaperPositionMode ='auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(['CRateGradient_',num2str(j)],'-dpdf');

% Polyfit the function of h and Gin_error to a first-order.
    p2 = polyfit(log10(h), log10(Gin_error), 1);
    % The convergence rate is the first number of 'p2'.
    convergenceRate2 = p2(1);
    % Display the value of convergence rate
    disp(['Convergence Rate of Gradient for interior nodes : ' num2str(convergenceRate2)]);
    % Calculate the convergence rate
    convergence_rate2 = [0
        diff(log(Gin_error)) ./ diff(log(h))];

    % Polyfit the function of h and Gh_error to a first-order.
    p3 = polyfit(log10(h), log10(Gh_error), 1);
    % The convergence rate is the first number of 'p3'.
    convergenceRate3 = p3(1);
    % Display the value of convergence rate
    disp(['Convergence Rate of Gradient : ' num2str(convergenceRate3)]);
    % Calculate the convergence rate
    convergence_rate3 = [0
        diff(log(Gh_error)) ./ diff(log(h))];


    % Polyfit the function of h and Hin_error to a first-order.
    p4 = polyfit(log10(h), log10(Hin_error), 1);
    % The convergence rate is the first number of 'p4'.
    convergenceRate4 = p4(1);
    % Display the value of convergence rate
    disp(['Convergence Rate of Hessian for interior nodes: ' num2str(convergenceRate4)]);
    % Calculate the convergence rate
    convergence_rate4 = [0
        diff(log(Hin_error)) ./ diff(log(h))];

    % Polyfit the function of h and Hh_error to a first-order.
    p5 = polyfit(log10(h), log10(Hh_error), 1);
    % The convergence rate is the first number of 'p5'.
    convergenceRate5 = p5(1);
    % Display the value of convergence rate
    disp(['Convergence Rate of Hessian : ' num2str(convergenceRate5)]);
    % Calculate the convergence rate
    convergence_rate5 = [0
        diff(log(Hh_error)) ./ diff(log(h))];

    % disp('Dofs  e_0  rate  e_1 rate e_2  rate e_3 rate e_4  rate');
    % for i = 1:length(h)
    %     fprintf('%.d & %.4e & %.2f & %.4e & %.2f & %.4e & %.2f & %.4e & %.2f & %.4e & %.2f  \n', ...
    %         Dofs(i), error(i), convergence_rate1(i),Gh_error(i), convergence_rate3(i),Gin_error(i), convergence_rate2(i), ...
    %         Hh_error(i), convergence_rate5(i),Hin_error(i), convergence_rate4(i));
    % end
        disp('Dofs  e_1 rate e_2  rate ');
    for i = 1:length(h)
        fprintf('%.d & %.4e & %.2f & %.4e & %.2f \n', ...
            Dofs(i),Gh_error(i), convergence_rate3(i),Gin_error(i), convergence_rate2(i));
    end

    disp('Dofs e_3 rate e_4  rate');
    for i = 1:length(h)
        fprintf('%.d & %.4e & %.2f & %.4e & %.2f  \n', ...
            Dofs(i), Hh_error(i), convergence_rate5(i),Hin_error(i), convergence_rate4(i));
    end
    disp('--------------------------------------------------------------------------------------');

end