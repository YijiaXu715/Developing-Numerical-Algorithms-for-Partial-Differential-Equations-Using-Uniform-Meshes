function [y,error]=rungekutta2(t0,tf,y0,h,ft_,fy_,m2)
% Backward Euler can solve a first-order ODE
% ft is the function of dy/dt=f(t,y)
% ft1 is the function of f'(t,y)
% Output y and error

% The number of time steps
% Rounding the substraction, abs=absolute value
N=round(abs(tf-t0)/h);

% Initialize x and y with zeros
t=zeros(N+1,1);
y=zeros(N+1,1);
t(1)=t0;
y(1)=y0;
p=1./(2*m2); 
m1=1-m2;

for i=1:N
    t(i+1)=t(i)+h;
    k1=ft_(t(i),y(i));
    k2=ft_(t(i)+p*h,y(i)+p*h*k1);
    y(i+1)=y(i)+h*(m1*k1+m2*k2);
end
yexact= fy_(t);

% The maximum of the error
error = max(abs(y-yexact));

% Numerical Solution
% figure(1)
% plot(t, y, '-o','Displayname','Numerical Solution');
% xlabel('t');
% ylabel('y');
% title('Numerical Solution by Second Order Runge Kutta Method');
% fig=gcf;
% fig.PaperPositionMode='auto';
% fig_pos=fig.PaperPosition;
% fig.PaperSize=[fig_pos(3) fig_pos(4)];
% print('NumericalSolution_1','-dpdf');

% Exact Solution
% figure(2)
% plot(t, y,'-o','Displayname','Numerical Solution');
% hold on
plot(t,yexact,'-*','Displayname','Exact Solution');
% hold off
% xlim([2.9,3.1]);
xlabel('t');
% ylim([2.4,2.8]);
ylabel('y');
% title('Exact Solution by computing');
legend('Location','best');
fig=gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];
print('comparison_2_h0','-dpdf');
end