function [u_h,Gux,Guy,Hxx,Hxy,Hyx,Hyy,u,ux,uy,uxx,uxy,uyx,uyy,inNode,node,elem] = elliptic_difference(nx,ny,degree)
[node,elem]=cartesian(nx,ny);
[Gux,Guy,Hxx,Hxy,Hyx,Hyy]=hessianrecovery2(node,elem);
Num_node = length(node);
x = node(:,1);
y = node(:,2);
[u,ux,uy,uxx,uxy,uyx,uyy] = exact_function(x,y,degree);
u_h = zeros(Num_node,1);
[bdEdge] = auxstructure_cartesian(elem);
bdNode = unique(bdEdge(:));
AllNode = 1:Num_node;
AllNode = AllNode';
inNode = setdiff(AllNode,bdNode);
%num_in = length(inNode);
f = -uxx-uyy;
% Evaluate the u_h of boundary nodes
u_h(bdNode) = u(bdNode);
% Evaluate the matrix f only with respect to the interior nodes
coefficient1 = -Hxx-Hyy;
f_rhs = f-coefficient1*u_h;
f_rhs = f_rhs(inNode(:));
% the coefficient of the interior nodes
coefficient = coefficient1(inNode(:),inNode(:));
uh = coefficient\f_rhs;
u_h(inNode) = uh; 
figure(1)
trisurf(elem,node(:,1),node(:,2),u_h,'EdgeColor','interp','FaceColor','none');
end