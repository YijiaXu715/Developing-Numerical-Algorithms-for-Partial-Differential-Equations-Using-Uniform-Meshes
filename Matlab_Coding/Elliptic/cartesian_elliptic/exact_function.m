function [u,ux,uy,uxx,uxy,uyx,uyy] =exact_function(x,y,degree)
N = size(x,1);

switch degree
    %% Smooth function
    case '1'
        u = sin(pi*x).*sin(pi*y);
        ux = pi*cos(pi*x).*sin(pi*y);
        uy = pi*sin(pi*x).*cos(pi*y);
        uxx = -pi^2*sin(pi*x).*sin(pi*y);
        uxy = pi^2*cos(pi*x).*cos(pi*y);
        uyx = pi^2*cos(pi*x).*cos(pi*y);
        uyy = -pi^2*sin(pi*x).*sin(pi*y);
        %% Quadratic
    case '2'
        u = x.^2+y.^2;
        ux = 2*x;
        uy = 2*y;
        uxx = 2*ones(N,1);
        uxy = zeros(N,1);
        uyx = zeros(N,1);
        uyy = 2*ones(N,1);
        %% Cubic
    case '3'
        u = 11*x.^3 + 37*x.^2.*y + 73*y.^2 + 87*x.*y + 71*ones(N,1);
        ux = 33*x.^2 + 74*x.*y + 87*y;
        uy = 37*x.^2 + 146*y + 87*x;
        uxx = 66*x + 74*y;
        uxy = 74*x + 87*ones(N,1);
        uyx = 74*x + 87*ones(N,1);
        uyy = 146*ones(N,1);
        %% Fourth
    case '4'
        u = x.^4 + 8*(x.^3).*y + 4*(x.^2).*(y.^2)+ 4*x.*(y.^3)+ 2*x.*y;
        ux = 4*(x.^3) + 24*(x.^2).*y + 8*x.*(y.^2)+ 4*(y.^3)+ 2*y;
        uy = 8*(x.^3) + 8*(x.^2).*y + 12*x.*(y.^2)+ 2*x;
        uxx = 12*x.^2 + 48*x.*y +8*y.^2;
        uyx = 24*x.^2 + 16*x.*y +12*y.^2 +2*ones(N,1);
        uxy = 24*x.^2 + 16*x.*y + 12*y.^2 +2*ones(N,1);
        uyy = 8*(x.^2) + 24*x.*y;
        %% Fifth
    case '5'
        u = 3*x.^5 + 4*y.^5 + 2*x.^3 + 3*y.^2;
        ux = 15*x.^4 + 6*x.^2;
        uy = 20*y.^4+6*y;
        uxx = 60*x.^3+12*x;
        uxy = zeros(N,1);
        uyx = zeros(N,1);
        uyy = 80*y.^3+6*ones(N,1);
        %% Sixth
    case '6'
        u = x.^6 + y.^6 + x.^2.*y.^4;
        ux = 6*x.^5 + 2*x.*y.^4;
        uy = 6*y.^5 + 4*x.^2.*y.^3;
        uxx = 30*x.^4 + 2*y.^4;
        uxy = 8*x.*y.^3;
        uyx = 8*x.*y.^3;
        uyy = 30*y.^4 + 12*x.^2.*y.^2;
        %% Seventh
    case '7'
        u = x.^7 + x.^5.*y.^2 + x.^3.*y.^3 + y.^7 ;
        ux = 7*x.^6 + 5*x.^4.*y.^2 + 3*x.^2.*y.^3 ;
        uy = 2*x.^5.*y + 3*x.^3.*y.^2 + 7*y.^6;
        uxx = 42*x.^5 + 20*x.^3.*y.^2 + 6*x.*y.^3;
        uxy = 10*x.^4.*y + 9*x.^2.*y.^2;
        uyx = 10*x.^4.*y + 9*x.^2.*y.^2;
        uyy = 2*x.^5 + 6*x.^3.*y + 42*y.^5;
end

