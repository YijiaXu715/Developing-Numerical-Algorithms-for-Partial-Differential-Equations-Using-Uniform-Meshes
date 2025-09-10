nx = 64;
ny = 64;
for j = 1:1
    if j==1
        disp(['Example' num2str(j) ':the example is a smooth function'])
    else
        disp(['Example ' num2str(j) ': the degree of the polynomial function is ' num2str(j)]);
    end
[u_h,Gux,Guy,Hxx,Hxy,Hyx,Hyy,u,ux,uy,uxx,uxy,uyx,uyy,inNode,node,elem] = elliptic_difference(nx,ny,num2str(j));
figure(1)
trisurf(elem,node(:,1),node(:,2),u_h,'EdgeColor','interp','FaceColor','none','FaceAlpha','0');
set(gca,'color','none'); 
fig=gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];
print(['Ru_h_',num2str(j)],'-dpng');

figure(2)
trisurf(elem,node(:,1),node(:,2),Gux*u_h,'EdgeColor','interp','FaceColor','none');
fig=gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];
print(['RGh_ux_',num2str(j)],'-dpdf');

figure(3)
trisurf(elem,node(:,1),node(:,2),Guy*u_h,'EdgeColor','interp','FaceColor','none');
fig=gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];
print(['RGh_uy_',num2str(j)],'-dpdf');

figure(4)
trisurf(elem,node(:,1),node(:,2),Hxx*u_h,'EdgeColor','interp','FaceColor','none');
fig=gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];
print(['RHh_xx_',num2str(j)],'-dpdf');

figure(5)
trisurf(elem,node(:,1),node(:,2),Hxy*u_h,'EdgeColor','interp','FaceColor','none');
fig=gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];
print(['RHh_xy_',num2str(j)],'-dpdf');

figure(6)
trisurf(elem,node(:,1),node(:,2),Hyy*u_h,'EdgeColor','interp','FaceColor','none');
fig=gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];
print(['RHh_yy_',num2str(j)],'-dpdf');
end