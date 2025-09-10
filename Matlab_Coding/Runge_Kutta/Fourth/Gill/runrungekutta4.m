% Calculate h for N times
N=4;
h = zeros(N,1);
error = zeros(N,1);
runtime = zeros(N,1);
for i=1:N
    % The cputime for iterations
    tic
    h(i) = 10^(-i);
    [y,error(i)]=rungekutta_gill(0,5,exp(-1),h(i),@ft,@fy);
    runtime(i) = toc;

    % Explanation of i,h,error and runtime
    disp(['i=', num2str(i), ', the step size is ' num2str(h(i)), ', the error is ' num2str(error(i))]);
    disp([ num2str(h(i)) ,'&', num2str(error(i))]);
    disp(['runtime is ', num2str(runtime(i))]);
end
figure(3);
loglog(1./h,error,'-o');
print('convergence_rate_n4','-dpdf');
xlabel('Numbers of Time steps');
ylabel('Errors');
title('Convergence Rate Analysis of Runge Kutta Gill Fourth order Euler Method');

% Polyfit the function of h and error to a first-order.
p = polyfit(log10(h), log10(error), 1);

% The convergence rate is the first number of 'p'.
convergenceRate = -p(1);

% Display the value of convergence rate
disp(['Convergence Rate of Runge Kutta Gill Fourth Order Method: ' num2str(convergenceRate)]);

% Calculate the convergence rate
convergence_rate = [0 
    diff(log(error)) ./ diff(log(h))];

disp('Step Size     Error        Convergence Rate      Runtime ');
disp('---------------------------------------------------');
for i = 1:length(h)
   fprintf('%.2e  & %.2e  & %.2f  & %.4f \n', h(i), error(i), convergence_rate(i),runtime(i));
end


