% Calculate h for N times
N=3;
t0=1;
h=0.001;
y0=sin(4)+1;
m2=0.75;
tf = zeros(N,1);
error=zeros(N,1);
runtime = zeros(N,1);


for i=1:N
    % The cputime for iterations
    tic
    tf(i) = 10^(i);
    [y,error(i)]=rungekutta2(t0,tf(i),y0,h,@ft_,@fy_,m2);
    runtime(i) = toc;

    % Explanation of i,h,error and runtime
    disp(['i=', num2str(i) ', the tf is ' num2str(tf(i))]);
    % disp([ num2str(h(i)) ,'&', num2str(error(i))]);
    %disp(['runtime is ', num2str(runtime(i))]);
end
figure(3);
loglog(tf,-1./error,'-o');
fig=gcf;
fig.PaperPositionMode='auto';
fig_pos=fig.PaperPosition;
fig.PaperSize=[fig_pos(3) fig_pos(4)];
print('runrktf','-dpdf');
xlabel('The final t (tf)');
ylabel('Error');
% title('Convergence Rate Analysis of Runge Kutta Euler Method');

% Polyfit the function of h and error to a first-order.
% p = polyfit(log10(tf), log10(error), 1);

% % The convergence rate is the first number of 'p'.
% convergenceRate = -p(1);
% 
% % Display the value of convergence rate
% disp(['Convergence Rate of Runge Kutta Method: ' num2str(convergenceRate)]);
% 
% % Calculate the convergence rate
% convergence_rate = [0 
%     diff(log(error)) ./ diff(log(h))];

disp('Final t   &  Error    &   Runtime ');
disp('---------------------------------------------------');
for i = 1:length(tf)
   fprintf('%.4e  & %.4e   & %.4f \n', tf(i), error(i),runtime(i));
end









