function [y,error]=rungekutta(t0,tf,y0,h,ft,fy,p,l)
% ft is the function of dy/dt=f(t,y)
% ft1 is the function of f'(t,y)
% fy is the exact solution 
% Output y and error

% The number of time steps
% Rounding the substraction, abs=absolute value
N=round(abs(tf-t0)/h);

% Initialize x and y with zeros
t=zeros(N+1,1);
y=zeros(N+1,1);
t(1)=t0;
y(1)=y0;
m2=(1/3-1/2*l)./(p.^2-p*l);
m3=(1/2*p-1/3)./(p*l-l.^2);
q=1./(6*p*m3);
r=l-q;
m1=1-m2-m3;

for i=1:N
    t(i+1)=t(i)+h;
    k1=ft(t(i),y(i));
    k2=ft(t(i)+p*h,y(i)+h*p*k1);
    k3=ft(t(i)+l*h,y(i)+h*r*k1+h*q*k2);
    y(i+1)=y(i)+h*(m1*k1+m2*k2+m3*k3);
end
yexact= fy(t);

% The maximum of the error
error = max(abs(y-yexact));

% Numerical Solution
% figure(1)
% plot(t, y, '-o','Displayname','Numerical Solution');
% xlabel('t');
% ylabel('y');
% title('Numerical Solution by Third Order Runge Kutta Method');
% fig=gcf;
% fig.PaperPositionMode='auto';
% fig_pos=fig.PaperPosition;
% fig.PaperSize=[fig_pos(3) fig_pos(4)];
% print('NumericalSolution_1','-dpdf');

% Exact Solution
figure(2)
plot(t, y,'-o','Displayname','Numerical Solution');
hold on
plot(t,yexact,'-*','Displayname','Exact Solution');
hold off
xlabel('t');
ylabel('y');
% title('Exact Solution by computing');
legend('Location','best');
fig=gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];
print('Comparison_NEh01','-dpdf');
end