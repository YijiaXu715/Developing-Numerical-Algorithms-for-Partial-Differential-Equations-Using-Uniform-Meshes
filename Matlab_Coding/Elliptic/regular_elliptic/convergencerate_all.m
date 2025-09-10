M = 5;
error = zeros(M,1);
Gh_error = zeros(M,1);
Gin_error = zeros(M,1);
Hh_error = zeros(M,1);
Hin_error = zeros(M,1);
h = zeros(M,1);
runtime = zeros(M,1);
Dofs = zeros(M,1);
for j = 1:7
    if j==1
        disp(['Example' num2str(j) ':the example is a smooth function'])
    else
        disp(['Example ' num2str(j) ': the degree of the polynomial function is ' num2str(j)]);
    end
    for i=1:M
        N = 2^(i+3);
        h(i) = 1/(2^(i+3));
        tic
        [u_h,Gux,Guy,Hxx,Hxy,Hyx,Hyy,u,ux,uy,uxx,uxy,uyx,uyy,inNode,node,elem] = elliptic_difference(N,N,num2str(j));
        Dofs(i) = size(u_h,1);
        % error for u
        error(i) = max(abs(u-u_h));
        % error for gradient recovery
        Gh_ux = Gux*u_h;
        Gh_uy = Guy*u_h;
        errorx =abs(ux-Gh_ux);
        errory= abs(uy-Gh_uy);
        L = 0.1;
        interior = find(node(:,1)<1-L&node(:,1)>L&node(:,2)<1-L&node(:,2)>L);
        in_errorX = max(errorx(interior(:,1),1));
        in_errorY = max(errory(interior(:,1),1));
        Gh_error(i) = max(max(errorx),max(errory));
        Gin_error(i) = max(in_errorX,in_errorY);
        % error for hr
        Hh_uxx = Hxx*u_h;
        Hh_uxy = Hxy*u_h;
        Hh_uyx = Hyx*u_h;
        Hh_uyy = Hyy*u_h;
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
    figure(1);
    loglog(Dofs,error,'-o');
    xlabel('Number of Degree of Freedoms');
    ylabel({'$\|u-u_h\|_{\infty,\Omega_h}$'},'Interpreter','latex');
    % fig = gcf;
    % fig.PaperPositionMode ='auto';
    % fig_pos = fig.PaperPosition;
    % fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(['Rate_u_',num2str(j)],'-dpdf');

    figure(2);
    loglog(Dofs,Gin_error,'-o');
    xlabel('Number of Degree of Freedoms');
    ylabel({'$\|G_hu_h-\nabla u\|_{\infty,\Omega_{1,h}}$'},'Interpreter','latex');
    % fig = gcf;
    % fig.PaperPositionMode ='auto';
    % fig_pos = fig.PaperPosition;
    % fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(['Rate_Gin_',num2str(j)],'-dpdf');

    figure(3);
    loglog(Dofs,Gh_error,'-o');
    xlabel('Number of Degree of Freedoms');
    ylabel({'$\|G_hu_h-\nabla u\|_{\infty,\Omega_h}$'},'Interpreter','latex');
    % fig = gcf;
    % fig.PaperPositionMode ='auto';
    % fig_pos = fig.PaperPosition;
    % fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(['Rate_Gh_',num2str(j)],'-dpdf');

    figure(4);
    loglog(Dofs,Hin_error,'-o');
    xlabel('Number of Degree of Freedoms');
    ylabel({'$\|H_hu_h-\Delta u\|_{\infty,\Omega_{1,h}}$'},'Interpreter','latex');
    % fig = gcf;
    % fig.PaperPositionMode ='auto';
    % fig_pos = fig.PaperPosition;
    % fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(['Rate_Hin_',num2str(j)],'-dpdf');

    figure(5);
    loglog(Dofs,Hh_error,'-o');
    xlabel('Number of Degree of Freedoms');
    ylabel({'$\|H_hu_h-\Delta u\|_{\infty,\Omega_h}$'},'Interpreter','latex');
    % fig = gcf;
    % fig.PaperPositionMode ='auto';
    % fig_pos = fig.PaperPosition;
    % fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(['Rate_Hh_',num2str(j)],'-dpdf');

    figure(6)
    loglog(Dofs,error,'k-o',Dofs,Gin_error,'r-*',Dofs,Gh_error,'b-*', ...
        Dofs,Hin_error,'m-+',Dofs,Hh_error,'-+');
    legend('$\|u-u_h\|_{\infty,\Omega_h}$','$\|G_hu_h-\nabla u\|_{\infty,\Omega_{1,h}}$','$\|G_hu_h-\nabla u\|_{\infty,\Omega_h}$', ...
        '$\|H_hu_h-\Delta u\|_{\infty,\Omega_{1,h}}$','$\|H_hu_h-\Delta u\|_{\infty,\Omega_h}$','Interpreter','latex','Location','best');
    xlabel('Number of Degree of Freedoms');
    ylabel('Error');
    % fig = gcf;
    % fig.PaperPositionMode ='auto';
    % fig_pos = fig.PaperPosition;
    % fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(['Rate_',num2str(j)],'-dpdf');

    % Polyfit the function of h and error to a first-order.
    p1 = polyfit(log10(h), log10(error), 1);
    % The convergence rate is the first number of 'p1'.
    convergenceRate1 = p1(1);
    % Display the value of convergence rate
    disp(['Convergence Rate of Error : ' num2str(convergenceRate1)]);
    % Calculate the convergence rate
    convergence_rate1 = [0
        diff(log(error)) ./ diff(log(h))];
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

    disp('Dofs  e_0  rate  e_1 rate e_2  rate e_3 rate e_4  rate');
    for i = 1:length(h)
        fprintf('%.d & %.2e & %.2f & %.2e & %.2f & %.2e & %.2f & %.2e & %.2f & %.2e & %.2f  \n', ...
            Dofs(i), error(i), convergence_rate1(i),Gh_error(i), convergence_rate3(i),Gin_error(i), convergence_rate2(i), ...
            Hh_error(i), convergence_rate5(i),Hin_error(i), convergence_rate4(i));
    end
    disp('--------------------------------------------------------------------------------------');
end

