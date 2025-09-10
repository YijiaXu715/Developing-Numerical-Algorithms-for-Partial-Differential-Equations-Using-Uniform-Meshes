function [x1]=functionnewton(f,g,x0,max_interations,epsilon)
for j=1:max_interations
    x1=x0-f(x0)./g(x0);
    if abs(f(x1))<epsilon 
        % The absolute value of the function is smaller than epsilon(error)
        % Satisfying the condition, it enters into the 'if' loop
       break
    else 
        x0=x1;
    end
end
end 

