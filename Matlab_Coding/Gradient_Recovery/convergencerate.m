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
    Gh_error = zeros(M,1);
    Gin_error = zeros(M,1);
    h = zeros(M,1);
    runtime = zeros(M,1);
    Dofs = zeros(M,1);
    for i=1:M
        N = 2^(i+3);
        h(i) = 1/(2^(i+3));
        tic
        [node,elem] = regular(N,N);
        Dofs(i) = size(node,1);
        [Gux,Guy] = gradientrecovery(node,elem);
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
        runtime(i) = toc;
    end
    figure(3);
    loglog(Dofs,Gh_error,'-o');
    xlabel('Number of Degree of Freedoms');
    ylabel('Error');
    fig = gcf;
    fig.PaperPositionMode ='auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(num2str(j),'-dpdf','-bestfit');

    % Polyfit the function of h and error to a first-order.
    p1 = polyfit(log10(h), log10(Gh_error), 1);
    % The convergence rate is the first number of 'p1'.
    convergenceRate1 = p1(1);
    % Display the value of convergence rate
    disp(['Convergence Rate : ' num2str(convergenceRate1)]);
    % Calculate the convergence rate
    convergence_rate1 = [0
        diff(log(Gh_error)) ./ diff(log(h))];

    disp('Dofs  Error   Convergence Rate of all   Runtime ');
    disp('---------------------------------------------------');
    for i = 1:length(h)
        fprintf('%.d & %.2e & %.2f & %.2f \n', Dofs(i), Gin_error(i), convergence_rate1(i),runtime(i));
    end


    figure(4);
    loglog(Dofs,Gin_error,'-o');
    xlabel('Number of Degree of Freedoms');
    ylabel('Error');
    fig = gcf;
    fig.PaperPositionMode ='auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(num2str(j),'-dpdf','-bestfit');

     % Polyfit the function of h and error to a first-order.
    p2 = polyfit(log10(h), log10(Gin_error), 1);
    % The convergence rate is the first number of 'p2'.
    convergenceRate2 = p2(1);
    % Display the value of convergence rate
    disp(['Convergence Rate : ' num2str(convergenceRate2)]);
    % Calculate the convergence rate
    convergence_rate2 = [0
        diff(log(Gin_error)) ./ diff(log(h))];

    disp('Dofs  Error   Convergence Rate of inner   Runtime ');
    disp('---------------------------------------------------');
    for i = 1:length(h)
        fprintf('%.d & %.2e & %.2f & %.2f \n', Dofs(i), Gh_error(i), convergence_rate2(i),runtime(i));
    end
end