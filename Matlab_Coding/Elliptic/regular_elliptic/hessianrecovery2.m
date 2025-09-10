function [Gux,Guy,Hxx,Hxy,Hyx,Hyy]=hessianrecovery2(node,elem)
Num_node=length(node);
Patch_elem=cell(Num_node,1);
% elem(i,1),elem(i,2),elem(i,3) let i plug into the proper row in cell
for i=1:length(elem)
    Patch_elem{elem(i,1)}(end+1)=i;
    Patch_elem{elem(i,2)}(end+1)=i;
    Patch_elem{elem(i,3)}(end+1)=i;
    % Patch_elem{elem(i,4)}(end+1)=i;
end
%% layer 1
Patch_node1 = cell(Num_node,1);
for i=1:Num_node
    % Go through Patch_elem{i} in the expression of elem
    % express elem in 'node'
    patch_node =elem((Patch_elem{i}(:)),:);
    % Transform the node elem in row matrix
    patch_node = patch_node(:)';
    % Cancel the same number in the node_elem
    patch_node =unique(patch_node);
    Patch_node1{i}=patch_node;
end
layer1=Patch_node1;
clear Patch_node1
%%  layer 2
Patch_node2 = cell(Num_node,1);
for i = 1:Num_node
    % Patch_node2{i}= patch_node;% the i-th row in Node_elem is node_elem
    if length(layer1{i}) > 6 % evaluate the number of the nodes in each rows
        Patch_node2{i} = layer1{i}; % if the size > 6, ends the process
    else % if the size < 6, extend to layer 2
        extendLayer2 = [];
        for j = 1:length(layer1{i}) % rotate 'j' times, 'j' is the number of columns of 'Patch_node'
            % plug 'nodes' layer 1 into Patch_node1
            extend_layer2 =layer1{layer1{i}(j)};
            % combine all layer 1 together
            extendLayer2(end+1:end+size(extend_layer2,2)) = extend_layer2;
            extendLayer2 = unique(extendLayer2);
        end
        Patch_node2{i} = extendLayer2;
    end
end
layer2 = Patch_node2;
clear Patch_node2
%% layer 3
Patch_node3 = cell(Num_node,1);
for i = 1:Num_node
    if length(layer1{i}) > 6 % evaluate the number of the nodes in each rows
        Patch_node3{i} = layer1{i}; % if the size > 6, ends the process
    else % if the size < 6, extend to layer 2
        extendLayer3 = [];
        if length(layer2{i}) > 8 % evaluate the number of the nodes in each rows
            Patch_node3{i} = layer2{i}; % if the size > 8, ends the process
        else % if the size < 8, extend to layer 3
            for j = 1:length(layer2{i})
                % plug 'nodes' of layer 2 after 'nodes' of layer 1
                extend_layer3 =layer1{layer2{i}(j)};
                % combine all layer 3
                extendLayer3(end+1:end+size(extend_layer3,2)) = extend_layer3;
                % cancel the same 'nodes'
                extendLayer3 = unique(extendLayer3);
            end
            % Combine the 'nodes' both centering and boundary
            Patch_node3{i} = extendLayer3;
        end
    end
end
layer3 = Patch_node3;
clear Patch_node3
%% Gradient Recovery
A1 = cell(Num_node,1);
indRow = [];
indCol = [];
gux = [];
guy = [];
hxx = [];
hxy = [];
hyx = [];
hyy = [];
for i=1:Num_node
    nodeCoordinate = node(layer3{i},:);
    max_x = max(abs(node(i,1)-nodeCoordinate(:,1)));
    max_y = max(abs(node(i,2)-nodeCoordinate(:,2)));
    h = max([max_x,max_y]);
    nodeCoordinate(:,1) = (nodeCoordinate(:,1)-node(i,1))/h;
    nodeCoordinate(:,2) = (nodeCoordinate(:,2)-node(i,2))/h;
    A = ones(length(nodeCoordinate),6);
    A(:,[2,3]) = nodeCoordinate;
    A(:,5) = nodeCoordinate(:,1).*nodeCoordinate(:,2);
    A(:,[4,6]) = nodeCoordinate.*nodeCoordinate;
    A1{i} = (A'*A)\A';
    indRow(end+1:end+size(layer3{i},2)) = i*ones(1,size(layer3{i},2));
    indCol(end+1:end+size(layer3{i},2)) = layer3{i};
    gux(end+1:end+size(layer3{i},2)) = A1{i}(2,:)/h;
    guy(end+1:end+size(layer3{i},2)) = A1{i}(3,:)/h;
    hxx(end+1:end+size(layer3{i},2)) = 2*A1{i}(4,:)/(h^2);
    hxy(end+1:end+size(layer3{i},2)) = A1{i}(5,:)/(h^2);
    hyx(end+1:end+size(layer3{i},2)) = A1{i}(5,:)/(h^2);
    hyy(end+1:end+size(layer3{i},2)) = 2*A1{i}(6,:)/(h^2);
end
Gux = sparse(indRow,indCol,gux);
Guy = sparse(indRow,indCol,guy);
Hxx = sparse(indRow,indCol,hxx);
Hxy = sparse(indRow,indCol,hxy);
Hyx = sparse(indRow,indCol,hyx);
Hyy = sparse(indRow,indCol,hyy);
end




