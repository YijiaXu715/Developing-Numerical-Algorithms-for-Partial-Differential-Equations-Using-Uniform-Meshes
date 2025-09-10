function [x1,j]=differencenewton(f,g,x0,max_interations,epsilon)
for j=1:max_interations
    x1=x0-f(x0)./g(x0);
    if abs(x1-x0)<epsilon
       break
    else 
        x0=x1;
    end
end
end      