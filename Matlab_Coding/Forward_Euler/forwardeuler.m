function [y,error] = forwardeuler(t0,tf,y0,h,rhs,exactSolution)
% Forward Euler can solve a first-order ODE
% rhs is the function of dy/dx=f(x,y)
% Output y and error

% The number of time steps
% Rounding the substraction, abs=absolute value
N = round(abs(tf-t0)/h);

% Initialize x and y with zeros
t = zeros(N+1,1);
y = zeros(N+1,1);
t(1) = t0;
y(1) = y0;

for i = 1:N
    t(i+1) = t(i)+h;
    y(i+1) = y(i)+h*rhs(t(i),y(i));
end
yexact = exactSolution(t);

% The maximum of the error
error = max(abs(y-yexact));

% Numerical Solution
figure(1)
plot(t, y, '-o','Displayname','Numerical Solution');
xlabel('x');
ylabel('y');
title('Numerical Solution by Forward Euler Method');
fig=gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];
print('NumericalSolution','-dpdf');

% % Exact Solution
% figure(2)
% plot(x, y, '-o','Displayname','Numerical Solution');
% hold on
% plot(x,yexact,'-*','Displayname','Exact Solution');
% hold off
% xlabel('x');
% ylabel('yexact');
% % title('Exact Solution by computing');
% legend('Location','best');
% fig=gcf;
% fig.PaperPositionMode='auto';
% fig_pos=fig.PaperPosition;
% fig.PaperSize=[fig_pos(3) fig_pos(4)];
% print('Solution_h2','-dpdf');

end