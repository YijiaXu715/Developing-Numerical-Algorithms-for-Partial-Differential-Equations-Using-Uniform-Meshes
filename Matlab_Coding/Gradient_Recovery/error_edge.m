clc
clear all
for j = 1:7
    if j==1
        disp(['Example' num2str(j) ':the example is a smooth function'])
    else
        disp(['Example ' num2str(j) ': the degree of the polynomial function is ' num2str(j)]);
    end
M = 4;
bd_error = zeros(M,1);
in_error = zeros(M,1);
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
    [u,ux,uy] = exactFunction(node,num2str(j));
    [bdEdge] = auxstructure_edge(elem);
    bdNode = unique(bdEdge(:));
    Gh_ux = Gux*u;
    Gh_uy = Guy*u;
   errorx = abs(ux-Gh_ux);
   errory = abs(uy-Gh_uy);
    %% Boundary
    bd_errorX = max(errorx(bdNode(:,1),1));
    bd_errorY = max(errory(bdNode(:,1),1));
    bd_error(i) = max(bd_errorX,bd_errorY);
    %% Interior
    Num_node = length(node);
    AllNode = 1:Num_node;
    AllNode = AllNode';
    inNode = setdiff(AllNode,bdNode);
    in_errorX = max(errorx(inNode(:,1),1));
    in_errorY = max(errory(inNode(:,1),1));
    in_error(i) = max(in_errorX,in_errorY);
    runtime(i) = toc;
end
disp('Dofs  bd_Error   in_Error   Runtime ');
disp('---------------------------------------------------');
for i = 1:length(h)
    fprintf('%.d & %.2e & %.2e & %.2f \n', Dofs(i), bd_error(i), in_error(i),runtime(i));
end
end
