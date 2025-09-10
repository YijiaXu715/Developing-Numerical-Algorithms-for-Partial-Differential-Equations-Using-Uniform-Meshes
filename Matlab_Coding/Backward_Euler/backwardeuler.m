function [y,error]=backwardeuler(t0,tf,y0,h,ft,fy,ft1,max_interations,epsilon)
% Backward Euler can solve a first-order ODE
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
    x0=y(i);
    f=@(x) -x+y(i)+h*ft(t(i+1),x); % Define the function of y(i+1)
    g=@(x) -1+h*ft1(t(i+1),x); % Define the derivative of y(i+1)
    % y(i+1)=functionnewton(f,g,x0,max_interations,epsilon);
    % y(i+1)=differencenewton(f,g,x0,max_interations,epsilon);
    y(i+1)=ratenewton(f,g,x0,max_interations,epsilon);
end
yexact= fy(t);

% The maximum of the error
error = max(abs(y-yexact));

% % Numerical Solution
% figure(1)
% plot(t, y, '-o','Displayname','Numerical Solution');
% xlabel('t');
% ylabel('y');
% title('Numerical Solution by Backward Euler Method');
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