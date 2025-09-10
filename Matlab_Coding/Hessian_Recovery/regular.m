function [node,elem]=regular(nx,ny)
% unit square
x0=0;
x1=1;
y0=0;
y1=1;
% the name of the function:regular
hx=(x1-x0)/nx;
hy=(y1-y0)/ny;
% hx,hy are the step of each partition
% nx,ny are the dimensions of the mesh
[x,y]=meshgrid(x0:hx:x1,y0:hy:y1);
% the matrix illustrates the point on the mesh
node=[x(:),y(:)];
% node: the coordinates of points on the mesh
elem=zeros(2*nx*ny,3);
% elem: show a null matrix with two columns
% plug the quantity of the elements and the quantity of the columns into it
for i=1:nx
    for j=1:ny
        elem(2*(j-1)*nx+2*(i-1)+1,:)=[j+(i-1)*(ny+1),j+1+i*(ny+1),j+1+(i-1)*(ny+1)];
        elem(2*(j-1)*nx+2*(i-1)+2,:)=[j+1+i*(ny+1),j+(i-1)*(ny+1),j+i*(ny+1)];
    end
end
% repeat the command by using different i,j
% elem(the number of elements,:)=[a,b,c]
% Connect three points that originate from the same point in the same direction,like anticlockwise
trisurf(elem,node(:,1),node(:,2),zeros(length(node),1),'facecolor','w');
view(2)
% the picture will show in two dimensions
fig=gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];
% adjust the papersize
print('regular','-dpdf','-r300');
% print (filename,filetype, resolution)
save('regular.mat','node','elem');
% save data
end