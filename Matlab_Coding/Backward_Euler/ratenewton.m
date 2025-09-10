function [x1]=ratenewton(f,g,x0,max_interations,epsilon)
for j=1:max_interations
    x1=x0-f(x0)./g(x0);
    if (abs(x1-x0))/x0<epsilon
        % The absolute value of the rate between x0 and x1 is smaller than epsilon(error)
        % Satisfying the condition, it enters into the 'if' loop
       break
    else 
        x0=x1;
    end
end
end      