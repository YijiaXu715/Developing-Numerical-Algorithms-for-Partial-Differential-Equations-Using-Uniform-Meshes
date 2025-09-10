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
    bd_error = zeros(M,1);
    in_error = zeros(M,1);
    h = zeros(M,1);
    runtime = zeros(M,1);
    Dofs = zeros(M,1);
    for i=1:M
        N = 2^(i+3);
        h(i) = 1/N;
        tic
        [node,elem] = regular(N,N);
        Dofs(i) = size(node,1);
        [Gux,Guy] = gradientrecovery(node,elem);
        [u,ux,uy] = exactFunction(node,num2str(j));
        Gh_ux = Gux*u;
        Gh_uy = Guy*u;
        errorx = abs(ux-Gh_ux);
        errory = abs(uy-Gh_uy);
         L = 0.1;
        interior = find(node(:,1)<1-L&node(:,1)>L&node(:,2)<1-L&node(:,2)>L);
        in_errorX = max(errorx(interior(:,1),1));
        in_errorY = max(errory(interior(:,1),1));
        in_error(i) = max(in_errorX,in_errorY);
        Num_node = length(node);
        AllNode = 1:Num_node;
        bdNode = (setdiff(AllNode,interior))';
        bd_errorX = max(errorx(bdNode(:,1),1));
        bd_errorY = max(errory(bdNode(:,1),1));
        bd_error(i) = max(bd_errorX,bd_errorY);
        runtime(i) = toc;
    end

    disp('Dofs  bd_error   in_error  Runtime ');
    disp('---------------------------------------------------');
    for i = 1:length(h)
        fprintf('%.d & %.2e & %.2e & %.2f \n', Dofs(i), bd_error(i), in_error(i),runtime(i));
    end
end