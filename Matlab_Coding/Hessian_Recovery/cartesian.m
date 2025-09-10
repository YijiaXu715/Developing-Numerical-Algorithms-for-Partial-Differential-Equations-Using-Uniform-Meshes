function [node,elem]=cartesian(nx,ny)
x0 = 0;
x1 = 1;
y0 = 0;
y1 = 1;
hx=(x1-x0)/nx;
hy=(y1-y0)/ny;
[x,y]=meshgrid(x0:hx:x1,y0:hy:y1);
node=[x(:),y(:)];
elem=zeros(nx*ny,4);
for i=1:nx
    for j=1:ny
        N=j+(i-1)*(ny+1);
        elem((j-1)*nx+i,:)=[N+ny+2,N+1,N,N+ny+1];
    end
end
% trisurf(elem,node(:,1),node(:,2),zeros(length(node),1),'facecolor','w');
% view(2);
% fig=gcf;
% fig.PaperPositionMode='auto';
% fig_pos=fig.PaperPosition;
% fig.PaperSize=[fig_pos(3) fig_pos(4)];
% print('cartesian','-dpdf','-r300');
% save('cartesian.mat','node','elem');
end