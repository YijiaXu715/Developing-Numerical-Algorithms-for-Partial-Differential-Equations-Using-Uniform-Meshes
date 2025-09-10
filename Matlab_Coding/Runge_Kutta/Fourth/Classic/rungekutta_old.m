function [y,error]=rungekutta_old(t0,tf,y0,h,ft,fy)
% ft is the function of dy/dt=f(t,y)
% ft1 is the function of f'(t,y)
% fy is the exact solution of f(t,y)
% Output y and error

% The number of time steps
% Rounding the substraction, abs=absolute value
N=round(abs(tf-t0)/h);

% Initialize x and y with zeros
t=zeros(N+1,1);
y=zeros(N+1,1);
t(1)=t0;
y(1)=y0;
for i=1:N
    t(i+1)=t(i)+h;
    k1=ft(t(i),y(i));
    k2=ft(t(i)+h./2,y(i)+h*k1./2);
    k3=ft(t(i)+h./2,y(i)+h*k2./2);
    k4=ft(t(i)+h,y(i)+h*k3);
    y(i+1)=y(i)+h*(k1+2*k2+2*k3+k4)./6;
end
yexact= fy(t);

% The maximum of the error
error = max(abs(y-yexact));

% % Numerical Solution
% figure(1)
% plot(t, y, '-o','Displayname','Numerical Solution');
% xlabel('t');
% ylabel('y');
% title('Numerical Solution by Fourth Order Runge Kutta Method');
% fig=gcf;
% fig.PaperPositionMode='auto';
% fig_pos=fig.PaperPosition;
% fig.PaperSize=[fig_pos(3) fig_pos(4)];
% print('NumericalSolution_1','-dpdf');
% 
% % Exact Solution
% figure(2)
% plot(t, y,'-o','Displayname','Numerical Solution');
% hold on
% plot(t,yexact,'-*','Displayname','Exact Solution');
% hold off
% xlabel('t');
% ylabel('yexact');
% title('Exact Solution by computing');
% legend('Location','best');
% fig=gcf;
% fig.PaperPositionMode='auto';
% fig_pos=fig.PaperPosition;
% fig.PaperSize=[fig_pos(3) fig_pos(4)];
% print('Comparison_NE1','-dpdf');
end